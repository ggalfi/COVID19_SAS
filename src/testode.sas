%ODEspec(
	dynvars=s1 i1 r1,
	dynfunc=%quote(
			deriv_s1=-s1*beta*i1;
			deriv_r1=i1*gamma;
			deriv_i1=-deriv_s1-deriv_r1;
			),
	outds=ode_rk4,
	method=RK4
);

%ODEspec(
	dynvars=s2 i2 r2,
	dynfunc=%quote(
			deriv_s2=-s2*beta*i2;
			deriv_r2=i2*gamma;
			deriv_i2=-deriv_s2-deriv_r2;
			),
	outds=ode_euler,
	method=euler
);


data testode;
	R0=3.0; gamma=0.05;	
	s1=0.99; i1=0.01; r1=0;
	s2=s1; i2=i1; r2=r1; 
	beta=R0*gamma;
	do day=1 to 100;
		%ODEstep(ode_rk4);
		%ODEstep(ode_euler);
		diffs=s1-s2; diffi=i1-i2; diffr=r1-r2; 
		output;		
	end;	
run;