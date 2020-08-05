

ods noproctitle;
ods graphics / imagemap=on;

proc sort data=COVID.MCMC_WEEKLY_POST out=WORK.TempSorted2236;
	by country region_code region_name;
run;

proc means data=WORK.TempSorted2236 chartype mean std min max median n 
		vardef=df qmethod=os;
	var mortality i0 R0_0 R0_1 beta gamma;
	output out=work.post_weekly_stat mean=std=min=max=median=n= / autoname;
	by country region_code region_name;
run;

proc datasets library=WORK noprint;
	delete TempSorted2236;
run;

proc sort data=COVID.MCMC_WEEKLY_QUAD_POST out=WORK.TempSorted2237;
	by country region_code region_name;
run;

proc means data=WORK.TempSorted2237 chartype mean std min max median n 
		vardef=df qmethod=os;
	var mortality i0 R0_0 R0_half R0_1 gamma;
	output out=work.post_weekly_stat mean=std=min=max=median=n= / autoname;
	by country region_code region_name;
run;

proc datasets library=WORK noprint;
	delete TempSorted2237;
run;


/*Prepare input data for model scoring with all sampled paramter set (weekly, quadratic R0) */
PROC SQL;
   CREATE TABLE COVID.MCMC_INPUT_WQ_POST_JOIN AS 
   SELECT t1.country, 
          t1.region_code, 
          t1.region_name, 
          t1.week, 
          t1.weekidx, 
          t1.nweeks,
          t1.valid_from, 
          t1.valid_to,
          t1.new_deceased,
          t1.region_population,
          t2.Iteration, 
          t2.mortality, 
          t2.i0, 
          t2.R0_0, 
          t2.R0_half, 
          t2.R0_1, 
          t2.gamma
      FROM COVID.MCMC_WEEKLY_INPUT t1, COVID.MCMC_WEEKLY_QUAD_POST t2
      WHERE (t1.country = t2.country AND t1.region_code = t2.region_code AND t1.region_name = t2.region_name)
      ORDER BY t1.country,
               t1.region_code,
               t1.region_name,
               t2.Iteration,
               t1.week;
QUIT;

/*Scoring the data*/
data mcmc_wq_post_full_preds(
		drop = arrpdec: arrinf: arrsus: arrR0: w a b s i r  _:
	);
	set covid.MCMC_INPUT_WQ_POST_JOIN;
	by country region_code region_name iteration;
	array arrpdec[40];
	array arrinf[40];
	array arrsus[40];
	array arrR0[40];
	retain arrpdec arrinf arrsus arrR0;
	if first.iteration then do;		
		a=2*R0_0+2*R0_1-4*R0_half;
		b=4*R0_half-3*R0_0-R0_1;
		i=i0;
		s=1-i;
		r=0;
		do w=1 to nweeks;
			tau=(w-1)/(nweeks-1);
			R0=a*tau*tau+b*tau+R0_0;
			arrR0[w]=R0;
			beta=R0*gamma;
			prevr=r;
			%ODEstep(covid.sir_weekly_spec);
			arrpdec[w]=(mortality*(r-prevr));
			arrinf[w]=i;
			arrsus[w]=s;
			*put d= i= s= mortality= dr=;
		end;
	end;
	tau=(weekidx-1)/(nweeks-1);
	
	pred_dec_weekly_ratio=arrpdec[weekidx];
	pred_infected_ratio=arrinf[weekidx];
	pred_susceptible_ratio=arrsus[weekidx];
	R0=arrR0[weekidx];
	pred_cases_ratio=1-pred_susceptible_ratio;
	pred_dec_weekly=pred_dec_weekly_ratio*region_population;
	rel_dev=abs(
		(pred_dec_weekly-new_deceased)/
		(pred_dec_weekly)
	);
run;

%calchpd(
	mcmc_wq_post_full_preds,
	mcmc_wq_post_preds_R0,
	R0,
	byfields=country region_code region_name week
);

/*comparing predicted and actual daily deaths (weekly, quadratic R0)*/
data mcmc_wq_pred_with_orig;
	set covid.mcmc_weekly_input;
	set covid.mcmc_weekly_quad_predsum;
run;

%macro plotMCMCWQpred();
	ods graphics on / imagemap=on /*reset=index*/;
	*ods html path="&coviddir/results/pred" (URL=NONE) file="pred_daily_deaths.html";
	proc sql;
		create table plot_mcmc_wq_pred as
			select distinct 
				country,
				region_name,
				region_code,
				rank_region_by_cases,
				lockdown,
				valid_from,
				valid_to,
				max(week) as max_week format date.,
				min(week) as min_week format date.
			from mcmc_wq_pred_with_orig
			group by country, region_name, region_code
			order by country, rank_region_by_cases;
	quit;

	data _null_;
		set plot_mcmc_wq_pred end=last;
		call symput('pt_reg'||strip(put(_N_, best.)), kstrip(region_name));
		call symput('pt_reg_code'||strip(put(_N_, best.)), kstrip(region_code));
		call symput('pt_ctry'||strip(put(_N_, best.)), kstrip(country));
		call symput('pt_lockdown'||strip(put(_N_, best.)), strip(put(lockdown, best.)));
		call symput('pt_vf'||strip(put(_N_, best.)), strip(put(valid_from, best.)));
		call symput('pt_vt'||strip(put(_N_, best.)), strip(put(valid_to, best.)));
		call symput('pt_start'||strip(put(_N_, best.)), strip(put(min_week, best.)));
		call symput('pt_end'||strip(put(_N_, best.)), strip(put(max_week, best.)));
		if last then 
			call symput('pt_cnt', strip(put(_N_, best.)));
	run;
	%do i=1 %to &pt_cnt;

	ods graphics on / imagemap;
	title "Bayesian predictions for weekly deaths (quadratic R0)";
	title2 "Country: &&pt_ctry&i, Region: &&pt_reg&i";
	PROC sgplot DATA =mcmc_wq_pred_with_orig(where =(region_name="&&pt_reg&i"));
		band x=week upper=hpdupper lower=hpdlower / legendlabel="HPD interval (95%)" fillattrs=(color=lightpink);
		needle x=week  y=new_deceased / lineattrs=(thickness=3px) legendlabel="Actual";
		series x=week  y=mean/lineattrs=GraphFit2 legendlabel="Predicted";		
	RUN;
	title "Bayesian predictions for quadratic R0";
	title2 "Country: &&pt_ctry&i, Region: &&pt_reg&i";
	PROC sgplot DATA =mcmc_wq_post_preds_R0(where =(region_name="&&pt_reg&i"));
		band x=week upper=hpdupper lower=hpdlower / legendlabel="HPD interval (95%)" fillattrs=(color=lightblue);			
		series x=week  y=mean/lineattrs=GraphFit legendlabel="R0 Predicted";		
	RUN;
	%end;
	*ods html close;
%mend;
%plotMCMCWQpred();