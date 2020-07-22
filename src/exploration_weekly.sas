

%macro plotWeeklyDeceased();
	ods graphics on / imagemap=on reset=index;
	*ods html path="&coviddir/results/exploration" (URL=NONE) file="explore.html";
	
	proc sql;
		create table plot_weekly_region as
			select distinct country, region_name, region_code, rank_region_by_cases, lockdown, valid_from, valid_to
			from covid.mcmc_weekly_abt
			where rank_region_by_cases<=5
			order by country, rank_region_by_cases;
	quit;

	data attrmap;
		length fillcolor $32 id $4 value $8;
		id='myid';
		value='Lockdown';
		fillcolor='lightgoldenrodyellow';
		output;
		value='';
		fillcolor='white';
		output;
	run;

	data _null_;
		set plot_weekly_region end=last;
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
	data plot_weekly_data;
		set covid.mcmc_weekly_abt(where =(region_name="&&pt_reg&i"));
		if(week>=lockdown) then block="Lockdown";
		else block="";
	run;
	ods graphics on / imagemap;
	title "Weekly deaths";
	title2 "Country: &&pt_ctry&i, Region: &&pt_reg&i";
	PROC sgplot DATA =plot_weekly_data DATTRMAP=attrmap;
		block x=week block=block /
			attrid=myid
 			lineattrs=(thickness=0px);;
		needle x=week  y=new_deceased / lineattrs=(thickness=3px) legendlabel='Deaths per day (actual)';
		*series x=date  y=int_dec_avg_dot/lineattrs=GraphFit2  legendlabel='Deaths per day (smoothed)';

		refline &&pt_vf&i / 
			name="vrref"
			axis=x legendlabel="Validity range" 
			lineattrs=(pattern=shortdash color=green);
		refline &&pt_vt&i / 
			axis=x  
			lineattrs=(pattern=shortdash color=green);
	RUN;
	%end;
	*ods html close;
%mend;
 %plotWeeklyDeceased();