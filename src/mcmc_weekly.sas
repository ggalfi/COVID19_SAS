data mcmc_weekly_outl;
	set covid.mcmc_weekly_abt;
	by  country region_code region_name;
	retain lag1week lag2week lag1nd lag2nd;
	if first.region_name then do;
		outlier=0;
		lag1week=week;
	end;
	else if(lag1week=lag2week or lag1nd<10 or lag2nd<10) then
		outlier=0;
	else do;
		multnd=exp((week-lag1week)*(log(lag1nd)-log(lag2nd))/(lag1week-lag2week));
		expnd=lag1nd*multnd;		
		outlier=
			abs(expnd-new_deceased)>max(3*sqrt(abs(expnd)), 3*sqrt(abs(new_deceased))) and
			(new_deceased<min(lag1nd/(3*multnd), lag1nd*multnd/3) 
				or new_deceased>max(3*lag1nd/multnd, 3*lag1nd*multnd));		
	end;
	if(not outlier) then do;
		lag2week=lag1week;
		lag1week=week;			
		lag2nd=lag1nd;
		lag1nd=max(new_deceased,0);
	end;
run;
proc sql;
	create table covid.mcmc_weekly_input as 
		select *, 
			min(week) as first_week format=date.,
			max(week) as last_week format=date.,
			(max(week) - min(week))/7 +1 as nweeks,
			(week - min(week))/7 +1 as weekidx
		from mcmc_weekly_outl
		where week>=valid_from and week<=valid_to and  rank_region_by_cases<=5
			and (not outlier)
			/*and region_name='Lombardia' and country='US'*/
		group by country, region_code, region_name;
quit;
data _null_;
	init=0.01;	
	call symput('mort_init', init);
	call symput('mort_mu', log(init));
	call symput('mort_var', log(10)*0.5);

	init=0.01;
	call symput('i0_init', init);
	call symput('i0_mu', log(init));
	call symput('i0_var', log(10)*1.5);

	init=3;
	call symput('R0_init', init);
	call symput('R0_mu', log(init));
	call symput('R0_var', log(10)*0.2);
	
	call symput('log_R0_lin_mu', 0);
	call symput('log_R0_lin_var', log(2)/20);
	
	call symput('log_R0_quad_mu', 0);
	call symput('log_R0_quad_var', log(2)/400);

	init=0.05;
	call symput('gamma_init', init);
	call symput('gamma_mu', log(init));
	call symput('gamma_var', log(10)*0.3);
run;

%ODEspec(
	dynvars=s i r,
	dynfunc=%quote(
			deriv_s=-s*beta*i;
			deriv_r=i*gamma;
			deriv_i=-deriv_s-deriv_r;
			),
	outds=covid.sir_weekly_spec,
	deltat=7,
	method=rk4
);

ods graphics on / imagemap=off reset=index;
ods html path="&coviddir/results/mcmc_weekly" (URL=NONE) file="mcmc.html";
title "MCMC analysis of weekly regional data (RK4 solver)";
proc mcmc
	data=covid.mcmc_weekly_input
	diag=all
	seed=20200601 nmc=1000000
	nbi=50000
	ntu=50000
	thin= 500
	propcov=ind
	mchistory= detailed
	outpost=covid.mcmc_weekly_post
	monitor=(mortality i0 R0_0 R0_1 R0_split /*log_R0_lin log_R0_quad*/ beta gamma)
	plots(smooth)=all;
	by country region_code region_name;
	
	preddist outpred=covid.mcmc_weekly_pred  saveparm;
	ods output
		predsumint=covid.mcmc_weekly_predsum
		postsumint=covid.mcmc_weekly_postsumint;
	parms
		log_mort &mort_mu
		log_i0 &i0_mu
		/*log_R0 &R0_mu*/
		log_R0_0 &R0_mu
		log_R0_1 &R0_mu
		R0_split 0.5
		/*log_R0_lin &log_R0_lin_mu
		log_R0_quad &log_R0_quad_mu*/
		log_gamma &gamma_mu;
	prior log_mort ~ normal(&mort_mu, sd=&mort_var);
	prior log_i0  ~ normal(&i0_mu, sd=&i0_var);
	*prior log_R0 ~ normal(&R0_mu, sd=&R0_var);
	prior log_R0_0 ~ normal(&R0_mu, sd=&R0_var);
	prior log_R0_1 ~ normal(&R0_mu, sd=&R0_var);
	prior R0_split ~ uniform(0.1, 0.9);
	*prior log_R0_lin~ normal(&log_R0_lin_mu, sd=&log_R0_lin_var);
	*prior log_R0_quad~ normal(&log_R0_quad_mu, sd=&log_R0_quad_var);
	prior log_gamma ~ normal(&gamma_mu, sd=&gamma_var);
	beginnodata;
		array arrpdec[52];
		mortality=exp(log_mort);
		i0=exp(log_i0);
		R0_0=exp(log_R0_0);
		R0_1=exp(log_R0_1);
		*R0_1=exp(log_R0+nweeks*log_R0_lin+nweeks*nweeks*log_R0_quad);
		gamma=exp(log_gamma);
		R0_split_idx=floor(nweeks*R0_split);
		*put Rzero= beta= gamma= i0=;

		i=i0;
		s=1-i;
		r=0;
		do w=1 to nweeks;			
			*R0=exp(log_R0+w*log_R0_lin+w*w*log_R0_quad);
			if(w<R0_split_idx) then R0=R0_0; else R0=R0_1;
			beta=R0*gamma;
			prevr=r;
			%ODEstep(covid.sir_weekly_spec);
			arrpdec[w]=(mortality*(r-prevr));
			*put d= i= s= mortality= dr=;
		end;
	endnodata;

	*prob_deceased=(mortality*dr);
	*model_deceased=region_population*prob_deceased;
	*put days= int_dec_avg_dot= prob_deceased= model_deceased=;
	model new_deceased ~ binomial(region_population, p=arrpdec[weekidx]);
quit;
ods html close;
