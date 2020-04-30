data covid.mcmc_input;
	set covid.mcmc_abt(where=(date>=valid_from and date<=valid_to and  rank_region_by_cases<=3  /*and region_name='Île-de-France'  and country='ITA'*/));
run;

data _null_;
	init=0.01;	
	call symput('mort_init', init);
	call symput('mort_mu', log(init));
	call symput('mort_var', log(10)*1.0);

	init=0.01;
	call symput('i0_init', init);
	call symput('i0_mu', log(init));
	call symput('i0_var', log(10)*1.5);

	init=3;
	call symput('R0_init', init);
	call symput('R0_mu', log(init));
	call symput('R0_var', log(10)*0.4);

	init=0.05;
	call symput('gamma_init', init);
	call symput('gamma_mu', log(init));
	call symput('gamma_var', log(10)*0.5);
run;

ods graphics on;
ods graphics on / imagemap=off;
title "MCMC analysis of regions";
proc mcmc
	data=covid.mcmc_input
	diag=all
	seed=20200424 nmc=1000000
	nbi=50000
	ntu=50000
	thin= 500
	propcov=ind
	mchistory= detail
	outpost=covid.mcmc_post
	monitor=(mortality i0 R0 beta gamma)
	plots(smooth)=all;
	by country region_code region_name;
	
	preddist outpred=covid.mcmc_pred  saveparm;
	ods output predsumint=covid.mcmc_predsum postsumint=covid.mcmc_postsumint;
	parms log_mort &mort_mu log_i0 &i0_mu log_R0 &R0_mu log_gamma &gamma_mu;
	prior log_mort ~ normal(&mort_mu, sd=&mort_var);
	prior log_i0  ~ normal(&i0_mu, sd=&i0_var);
	
	prior log_R0 ~ normal(&R0_mu, sd=&R0_var);
	prior log_gamma ~ normal(&gamma_mu, sd=&gamma_var);
	beginnodata;
		days=valid_to-valid_from+1;
		array arrpdec[120];
		mortality=exp(log_mort);
		i0=exp(log_i0);
		R0=exp(log_R0);
		gamma=exp(log_gamma);
		beta=R0*gamma;
		*put Rzero= beta= gamma= i0=;

		i=i0;
		s=1-i;	
		do d=1 to days;
			ds=-s*beta*i;
			dr=i*gamma;
			s=s+ds;
			i=i-ds-dr;
			arrpdec[d]=(mortality*dr);
			*put d= i= s= mortality= dr=;
		end;
	endnodata;

	*prob_deceased=(mortality*dr);
	*model_deceased=region_population*prob_deceased;
	*put days= int_dec_avg_dot= prob_deceased= model_deceased=;
	model int_dec_avg_dot ~ binomial(region_population, p=arrpdec[date-valid_from+1]);
quit;

