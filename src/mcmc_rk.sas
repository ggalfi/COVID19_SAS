/*An alternative for mcmc.sas with a more accurate Runge-Kutta solver*/
%ODEspec(
	dynvars=s i r,
	dynfunc=%quote(
			deriv_s=-s*beta*i;
			deriv_r=i*gamma;
			deriv_i=-deriv_s-deriv_r;
			),
	outds=covid.sir_rk_spec,
	method=rk4
);

ods graphics on / imagemap=off;
ods html path="&coviddir/results/mcmc" (URL=NONE) file="mcmc_rk.html";
title "MCMC analysis of regions (Runge-Kutta solver)";
proc mcmc
	data=covid.mcmc_input
	diag=all
	seed=20200424 nmc=1000000
	nbi=50000
	ntu=50000
	thin= 500
	propcov=ind
	mchistory= detailed
	outpost=covid.mcmc_rk_post
	monitor=(mortality i0 R0 beta gamma)
	plots(smooth)=all;
	by country region_code region_name;
	
	preddist outpred=covid.mcmc_rk_pred  saveparm;
	ods output
		predsumint=covid.mcmc_rk_predsum
		postsumint=covid.mcmc_rk_postsumint;
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
		r=0;
		do d=1 to days;
			prevr=r;
			%ODEstep(covid.sir_rk_spec);
			arrpdec[d]=(mortality*(r-prevr));
		end;
	endnodata;

	model int_dec_avg_dot ~ binomial(region_population, p=arrpdec[date-valid_from+1]);
quit;
ods html close;