options mprint;
%macro ODEspec(dynvars, dynfunc, outds, deltat=1.0, timevar=, method=euler);
	%put &dynvars, &dynfunc;
	%local n;
	data _ode_tmp;
		length type 8 nval 8 cval $16384;
		type=1;
		nval=1;
		do while (scan("&dynvars", nval, ' ') ne '');
			cval=scan("&dynvars", nval, ' ');
			output;
			nval=nval+1;
		end;		
		type=0;
		cval='';
		nval=nval-1;
		output;
		call symput('n', strip(put(nval, best.)));
		type=2;
		nval=1;
		do while (scan("&dynfunc", nval, ';') ne '');
			cval=scan("&dynfunc", nval, ';');
			output;
			nval=nval+1;
		end;		
	run;
	data &outds;
		set _ode_tmp;
		length part $8 idx 8 line $16384;
		array varname_{&n} $32;
		retain varname_;
		n=&n;
		if(type=1);
		varname_[nval]=cval;
		if (nval=&n) then output;
		keep part idx line varname_:;
	run;
	%if (%upcase(&method) = EULER) %then %do;
	data &outds;
		set _ode_tmp(where=(type=2) in=intmp) end=last;
		if(_N_=1) then set &outds;
		array varname {*} varname_:;
		retain varname;
		retain  idx 0;
		if(intmp) then do;
			part='STEP';
			line=cval;
			idx+1;
			output;
		end;
		if last then do;
			part='STEP';
			do i=1 to &n;
				line=strip(varname[i])||' = '||strip(varname[i])||'+'||"(&deltat)*deriv_"||strip(varname[i]);
				idx+1;
				output;
			end;
		end;
		keep part idx line;
	run;
	%end;
	%else %if (%upcase(&method) = RK4) %then %do;
	data &outds;
		set
			_ode_tmp(where=(type=2) in=intmp1)
			_ode_tmp(where=(type=2) in=intmp2)
			_ode_tmp(where=(type=2) in=intmp3)
			_ode_tmp(where=(type=2) in=intmp4)
			end=last;
		retain prevpart 0;
		
		retain idx;
		
		if(_N_=1) then do;
			set &outds;
			array varname {*} varname_:;
			retain varname;
			idx=0;
			part='STEP';
			do i=1 to &n;
				line='_odetmp_0_'||strip(i)||'='||strip(varname[i]);
				idx+1;
				output;
			end;
		end;
		
		/*Calc k1*/
		if(intmp1) then do;
			part='STEP';
			line=cval;
			idx+1;
			output;
			prevpart=1;
		end;
		
		/*Calc k2*/
		if(intmp2) then do;
			if(prevpart=1) then do;
				do i=1 to &n;
					/*k1 saved to  _odetmp_1_<i>*/
					line='_odetmp_1_'||strip(i)||'= deriv_'||strip(varname[i]);
					idx+1;
					output;
			
					/*y+h*k1/2*/
					line=strip(varname[i])||' = '||strip(varname[i])||'+'||"0.5*(&deltat)*deriv_"||strip(varname[i]);
					idx+1;
					output;
				end;
			end;
			part='STEP';
			line=cval;
			idx+1;
			output;
			prevpart=2;
		end;
		
		/*Calc k3*/
		if(intmp3) then do;
			if(prevpart=2) then do;
				do i=1 to &n;
					/*k2 saved to  _odetmp_2_<i>*/
					line='_odetmp_2_'||strip(i)||'= deriv_'||strip(varname[i]);
					idx+1;
					output;
			
					/*y+h*k2/2*/
					line=strip(varname[i])||' = _odetmp_0_'||strip(i)||'+'||"0.5*(&deltat)*deriv_"||strip(varname[i]);
					idx+1;
					output;
				end;
			end;
			part='STEP';
			line=cval;
			idx+1;
			output;
			prevpart=3;
		end;
		
		/*Calc k4*/
		if(intmp4) then do;
			if(prevpart=3) then do;
				do i=1 to &n;
					/*k3 saved to  _odetmp_3_<i>*/
					line='_odetmp_3_'||strip(i)||'= deriv_'||strip(varname[i]);
					idx+1;
					output;
			
					/*y+h*k3*/
					line=strip(varname[i])||' = _odetmp_0_'||strip(i)||'+'||"(&deltat)*deriv_"||strip(varname[i]);
					idx+1;
					output;
				end;
			end;
			part='STEP';
			line=cval;
			idx+1;
			output;
			prevpart=4;
		end;
		
		if last then do;
			part='STEP';
			do i=1 to &n;
				line=strip(varname[i])||' = _odetmp_0_'||strip(i)||'+'||"((&deltat)/6)*("||
					' _odetmp_1_'||strip(i)||'+'||
					' 2*_odetmp_2_'||strip(i)||'+'||
					' 2*_odetmp_3_'||strip(i)||'+'||
					' deriv_'||strip(varname[i])||')'					
					;
				idx+1;
				output;
			end;
		end;
		keep part idx line;
	run;
	%end;
	%else %do;
		%put  No such ode method: &odetype;
		%abort;
	%end;
	/*proc sort data=_ode_tmp; by type nval;run;*/
%mend;

%macro ODEinternOpen(ds);
	%let _ode_dsid=%sysfunc(open(&dsspec,i));
	%if (&_ode_dsid = 0) %then %do;
   		%put %sysfunc(sysmsg());
   		%abort;
   	%end;
%mend;

%macro ODEinternClose(ds);
	%local rc;
   	%let rc=%sysfunc(close(&_ode_dsid));
%mend;

%macro ODEinternFindvar(varname);
	%local varid;
	%let varid=%sysfunc(varnum(&_ode_dsid, &varname));
   	%if ( &varid = 0 ) %then %do;
   		%put Variable &varname does not exist in dataset;
   		%sysfunc(close(&_ode_dsid));
   		%abort;
   	%end;
   	&varid
%mend;

%macro ODEstep(dsspec);
	%local _ode_dsid part line;
	
   	%ODEinternOpen(dsspec);
   	
   	%let part=%ODEinternFindvar(part);
   	%let line=%ODEinternFindvar(line);
   	%put _ode_dsid =&_ode_dsid;

   	%do %WHILE (%sysfunc(fetch(&_ode_dsid)) = 0);
   		%if( %sysfunc(getvarc(&_ode_dsid, &part)) = STEP ) %then %do;
			%sysfunc(getvarc(&_ode_dsid, &line));
		%end;
   	%end;
	
   	%ODEinternClose(dsspec);
%mend;


