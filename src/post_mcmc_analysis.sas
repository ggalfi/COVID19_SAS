/*Prepare input data for model scoring with all sampled paramter set */
PROC SQL;
   CREATE TABLE COVID.MCMC_INPUT_POST_JOIN AS 
   SELECT t1.country, 
          t1.region_code, 
          t1.region_name, 
          t1.date, 
          t1.valid_from, 
          t1.valid_to,
          t1.deceased_avg_dot,
          t1.region_population,
          t2.Iteration, 
          t2.mortality, 
          t2.i0, 
          t2.R0, 
          t2.beta, 
          t2.gamma
      FROM COVID.MCMC_INPUT t1, COVID.MCMC_POST t2
      WHERE (t1.country = t2.country AND t1.region_code = t2.region_code AND t1.region_name = t2.region_name)
      ORDER BY t1.country,
               t1.region_code,
               t1.region_name,
               t2.Iteration,
               t1.date;
QUIT;

/*Scoring the data*/
data 
	mcmc_post_full_preds(
		keep = country region_code region_name 
			iteration date pred_dec_daily_ratio
			pred_infected_ratio pred_susceptible_ratio
			pred_cases_ratio rel_dev
	);
	set covid.MCMC_INPUT_POST_JOIN;
	by country region_code region_name iteration;
	array arrpdec[120];
	array arrinf[120];
	array arrsus[120];
	array arrreldev[120];
	retain arrpdec arrinf arrsus;
	if first.iteration then do;
		days=valid_to-valid_from+1;	
		i=i0;
		s=1-i;	
		do d=1 to days;
			ds=-s*beta*i;
			dr=i*gamma;
			s=s+ds;
			i=i-ds-dr;
			arrpdec[d]=(mortality*dr);
			arrinf[d]=i;
			arrsus[d]=s;
		end;
	end;
	d=date-valid_from+1;
	pred_dec_daily_ratio=arrpdec[d];
	pred_infected_ratio=arrinf[d];
	pred_susceptible_ratio=arrsus[d];
	pred_cases_ratio=1-pred_susceptible_ratio;
	pred_dec_daily=pred_dec_daily_ratio*region_population;
	rel_dev=abs(
		(pred_dec_daily-deceased_avg_dot)/
		(pred_dec_daily)
	);
run;

proc sql;
	create table mcmc_post_rel_dev as
		select country, region_code, region_name, iteration,
 			mean(rel_dev) as rel_dev_mean, median(rel_dev) as rel_dev_median
		from mcmc_post_full_preds
		group by country, region_code, region_name, iteration;
quit;

%calchpd(
	mcmc_post_rel_dev,
	mcmc_rel_dev_stat_median,
	rel_dev_median,
	byfields=country region_code region_name
);

%calchpd(
	mcmc_post_rel_dev,
	mcmc_rel_dev_stat_mean,
	rel_dev_mean,
	byfields=country region_code region_name
);

/*Posterior statistics*/
proc sql;
	create table COVID.MCMC_POSTSUMINT_EXT as
		select * from (
			select * from COVID.MCMC_POSTSUMINT
			outer union corr
			select * from mcmc_rel_dev_stat_median 
			outer union corr
			select * from mcmc_rel_dev_stat_mean
		)
		order by country, region_code, region_name;			
run;

proc transpose data=COVID.MCMC_POSTSUMINT_EXT out=MCMC_POSTSUMINT_MEANS;
	by country region_code region_name;
	var mean;
	id parameter;
run;

proc sql noprint;
	select "'"||strip(region_name)||"'"
	into :selected_regions separated by ', '
	from MCMC_POSTSUMINT_MEANS
	where rel_dev_median<0.15;
run;
%put &selected_regions;

ods graphics on / imagemap=on reset=index;
ods html path="&coviddir/results/post_params" (URL=NONE) file="pred_params.html";
title "SIR parameters of regions";

proc sgplot data=MCMC_POSTSUMINT_MEANS;
	scatter x=rel_dev_mean y=rel_dev_median/ datalabel=region_name group=country;
	refline 0.14 / axis=y;
quit;

proc sgplot data=MCMC_POSTSUMINT_MEANS;
	scatter x=R0 y=mortality/ datalabel=region_name group=country;
quit;

proc sgplot data=MCMC_POSTSUMINT_MEANS;
	scatter x=R0 y=gamma/ datalabel=region_name group=country;
quit;

proc sgplot data=MCMC_POSTSUMINT_MEANS;
	scatter x=beta y=gamma/ datalabel=region_name group=country;
quit;

proc sgplot data=MCMC_POSTSUMINT_MEANS;
	scatter x=mortality y=gamma/ datalabel=region_name group=country;
quit;

proc sgplot data=MCMC_POSTSUMINT_MEANS;
	scatter x=rel_dev_median y=mortality/ datalabel=region_name group=country;
quit;
ods html close;

/*Overall spread of SIR parameters*/
%calchpd(
	covid.mcmc_post(where=(region_name in (&selected_regions))),
	mcmc_selected_stat_R0,
	R0
);

%calchpd(
	covid.mcmc_post(where=(region_name in (&selected_regions))),
	mcmc_selected_stat_mortality,
	mortality
);

data mcmc_selected_stat_mortality;
	set mcmc_selected_stat_mortality;
	format mean percent8.2;
	format stddev percent8.2;
	format hpdupper percent8.2;
	format hpdlower percent8.2;
run;

%calchpd(
	covid.mcmc_post(where=(region_name in (&selected_regions))),
	mcmc_selected_stat_gamma,
	gamma
);

%calchpd(
	covid.mcmc_post(where=(region_name in (&selected_regions))),
	mcmc_selected_stat_i0,
	i0
);



/*Posterior distribution of SIR paramteres*/
title "Prior and posterior distribution of R0";
PROC UNIVARIATE DATA = covid.mcmc_post
		CIBASIC(TYPE=TWOSIDED ALPHA=0.05)
		MU0=0
;
	by country;
	class  region_name ;
	VAR R0;
	HISTOGRAM R0 / nrows=3 nobars
		LOGNORMAL(W=1 L=1 COLOR=LIME ZETA=&R0_mu THETA=0 SIGMA=&R0_var)
		KERNEL(W=1 L=1 COLOR=CX008080 C=MISE K=NORMAL);
	ods select histogram;

RUN;

title "Prior and posterior distribution of mortality rate";
PROC UNIVARIATE DATA = covid.mcmc_post
		CIBASIC(TYPE=TWOSIDED ALPHA=0.05)
		MU0=0
;
	BY country;
	class region_name;
	VAR mortality;
	HISTOGRAM mortality/ nrows=3 nobars
		LOGNORMAL(W=1 L=1 COLOR=LIME ZETA=&R0_mu THETA=0 SIGMA=&R0_var)
		KERNEL(W=1 L=1 COLOR=CX008080 C=MISE K=NORMAL);
	ods select histogram;

RUN;

title "Prior and posterior distribution of recovery rate";
PROC UNIVARIATE DATA = covid.mcmc_post
		CIBASIC(TYPE=TWOSIDED ALPHA=0.05)
		MU0=0
;
	BY country;
	class region_name;
	VAR gamma;
	HISTOGRAM gamma/ nrows=3 nobars
		LOGNORMAL(W=1 L=1 COLOR=LIME ZETA=&R0_mu THETA=0 SIGMA=&R0_var)
		KERNEL(W=1 L=1 COLOR=CX008080 C=MISE K=NORMAL);
	ods select histogram;

RUN;
proc sql;
	create table mcmc_pred_uncertain as
		select
			country, 
			region_code, 
			region_name,
			mean((HPDUpper-HPDLower)/mean) as pred_unc_ratio
		from covid.mcmc_predsum
		group by country, region_code, region_name
		order by pred_unc_ratio;
quit;

/*comparing predicted and actual daily deaths */
data mcmc_pred_with_orig;
	set covid.mcmc_input;
	set covid.mcmc_predsum;
run;

%macro plotMCMCpred();
	ods graphics on / imagemap=on reset=index;
	ods html path="&coviddir/results/pred" (URL=NONE) file="pred_daily_deaths.html";
	title "Prediction for the daily deathes";
	proc sql;
		create table plot_mcmc_pred as
			select distinct country, region_name, region_code, rank_region_by_cases, lockdown, valid_from, valid_to
			from mcmc_pred_with_orig
			order by country, rank_region_by_cases;
	quit;

	data _null_;
		set plot_mcmc_pred end=last;
		call symput('pt_reg'||strip(put(_N_, best.)), kstrip(region_name));
		call symput('pt_reg_code'||strip(put(_N_, best.)), kstrip(region_code));
		call symput('pt_ctry'||strip(put(_N_, best.)), kstrip(country));
		call symput('pt_lockdown'||strip(put(_N_, best.)), strip(put(lockdown, best.)));
		call symput('pt_vf'||strip(put(_N_, best.)), strip(put(valid_from, best.)));
		call symput('pt_vt'||strip(put(_N_, best.)), strip(put(valid_to, best.)));
		if last then 
			call symput('pt_cnt', strip(put(_N_, best.)));
	run;
	%do i=1 %to &pt_cnt;

	ods graphics on / imagemap;
	title "Bayesian predictions for death per day";
	title2 "Country: &&pt_ctry&i, Region: &&pt_reg&i";
	PROC sgplot DATA =mcmc_pred_with_orig(where =(region_name="&&pt_reg&i"));
		band x=date upper=hpdupper lower=hpdlower / legendlabel="HPD interval (95%)" fillattrs=(color=lightpink);
		needle x=date  y=deceased_dot / lineattrs=(thickness=3px) legendlabel="Actual";
		series x=date  y=mean/lineattrs=GraphFit2 legendlabel="Predicted";		
	RUN;
	%end;
	ods html close;
%mend;
%plotMCMCpred();



/*Comparing the predictions with the reported cases*/
%calchpd(
	mcmc_post_full_preds,
	mcmc_post_preds_c,
	pred_cases_ratio,
	byfields=country region_code region_name date
);

PROC SQL;
   CREATE TABLE WORK.MCMC_INPUT_PREDS AS 
   SELECT t1.*, 
          t2.Mean as pred_cases_ratio_mean, 
          t2.HPDUpper as pred_cases_ratio_hpdu, 
          t2.HPDLower as pred_cases_ratio_hpdl
      FROM COVID.MCMC_INPUT t1, WORK.mcmc_post_preds_c t2
      WHERE (t1.country = t2.country AND t1.region_code = t2.region_code AND t1.region_name = t2.region_name AND t1.date 
           = t2.date)
      ORDER BY t1.country,
               t1.region_code,
               t1.region_name,
               t1.date;
QUIT;

data MCMC_INPUT_PREDS2;
	set MCMC_INPUT_PREDS;
	*pred_recovered_ratio=pred_susceptible_ratio+pred_infected_ratio;
	*pred_susceptible=pred_susceptible_ratio*region_population;
	*pred_infected=pred_infected_ratio*region_population;
	pred_total_cases=region_population*pred_cases_ratio_mean;
	pred_total_cases_hpdu=region_population*pred_cases_ratio_hpdu;
	pred_total_cases_hpdl=region_population*pred_cases_ratio_hpdl;
	format act_to_pred_total_cases percent8.3;
	act_to_pred_total_cases=total_cases/pred_total_cases;
run;

%macro plotMCMCtotcases();
	ods graphics on / imagemap=off reset=index(100);
	ods html path="&coviddir/results/pred" (URL=NONE) file="pred_total_cases.html";
	title "Prediction for the number of total cases";
	proc sql;
		create table plot_mcmc_pred_2 as
			select distinct country, region_name, region_code, rank_region_by_cases, lockdown, valid_from, valid_to
			from mcmc_pred_with_orig
			order by country, rank_region_by_cases;
	quit;

	data _null_;
		set plot_mcmc_pred_2 end=last;
		call symput('pt_reg'||strip(put(_N_, best.)), kstrip(region_name));
		call symput('pt_reg_code'||strip(put(_N_, best.)), kstrip(region_code));
		call symput('pt_ctry'||strip(put(_N_, best.)), kstrip(country));
		call symput('pt_lockdown'||strip(put(_N_, best.)), strip(put(lockdown, best.)));
		call symput('pt_vf'||strip(put(_N_, best.)), strip(put(valid_from, best.)));
		call symput('pt_vt'||strip(put(_N_, best.)), strip(put(valid_to, best.)));
		if last then 
			call symput('pt_cnt', strip(put(_N_, best.)));
	run;
	%do i=1 %to &pt_cnt;

	ods graphics on / imagemap;
	title "Total confirmed cases (Actual vs. Predicted)";
	title2 "Country: &&pt_ctry&i, Region: &&pt_reg&i";
	PROC sgplot DATA =MCMC_INPUT_PREDS2(where =(region_name="&&pt_reg&i"));
		band x=date upper=pred_total_cases_hpdu lower=pred_total_cases_hpdl / legendlabel="HPD interval (95%)" fillattrs=(color=lightgreen);
		series x=date  y=total_cases/y2axis lineattrs=GraphFit legendlabel="Actual";		
		series x=date  y=pred_total_cases/lineattrs=GraphFit2  lineattrs=(color=green) legendlabel="Predicted";		
	RUN;
	%end;
	ods html close;
%mend;
%plotMCMCtotcases();
