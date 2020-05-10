%macro calcHPD(inds, outds, var, byfields=, alpha=0.05);
	data _null_;
		length bycom $1024 byjoin $4096 lastfld $256;
		res="&byfields";
		i=1;
		do while(scan(res, i, ' ') ne '');
			if(i>1) then bycom=strip(bycom)||',';
			lastfld=scan(res, i, ' ');
			bycom=strip(bycom)||'h.'||lastfld;
			byjoin=strip(byjoin)||' and h.'||strip(lastfld)||' =  l.'||lastfld;
			i=i+1;
		end;
		call symput('bycom', strip(bycom));
		call symput('byjoin', strip(byjoin));
		call symput('lastfld', strip(lastfld));
	run;
	proc sql;
		create table _hpdtmp1 as
			select
	%if(&byfields ne ) %then &bycom,;
				&var,
				count(&var) as _hpdnmiss,
				mean(&var) as _hpdmean,
				std(&var) as _hpdstd
			from &inds as h
			where &var ne .
    %if(&byfields ne ) %then %do;
    		group by &bycom
    		order by &bycom
    %end;
    		;
	quit;
	
	proc rank data=_hpdtmp1 out=_hpdtmp1 ties=low;
	%if(&byfields ne ) %then %do;
		by &byfields;	
	%end;
		var &var;
		ranks _hpdidx;
	quit;
	
	data _hpdtmp1;
		set _hpdtmp1;
		_hpddelta=floor(_hpdnmiss*(1-&alpha));
		_hpdhiidx=_hpdidx+_hpddelta;
		if(_hpdhiidx<=_hpdnmiss) or (_hpdidx-_hpddelta>=1);
	run;
	
	proc sql;
		create table &outds as
			select
	%if(&byfields ne ) %then &bycom,;
				"&var" as Parameter length 32,
				h._hpdnmiss as N format 8.,
				h._hpdmean as Mean format d8.,
				h._hpdstd as StdDev format d8.,
				l.&var as HPDLower format d8.,
				h.&var as HPDUpper format d8.,
				(h.&var - l.&var) as HPDDiff
			from
				_hpdtmp1 as h
				inner join _hpdtmp1 as l
				on(h._hpdidx=l._hpdhiidx &byjoin);
	quit;
	proc sort data=&outds; by &byfields HPDDiff;run;
	
	
	data &outds;
	%if(&byfields ne ) %then %do;
		set &outds;
		by &byfields;
		if first.&lastfld;
	%end;
	%else %do;
		set &outds(obs=1);
	%end;
		drop HPDDiff;
	run;
	
	proc sql noprint;
		drop table _hpdtmp1;
	quit;
%mend;
*%calcHPD(COVID.MCMC_POST, testhpd, R0, byfields=country region_code region_name);
*%calcHPD(COVID.MCMC_POST, testhpd, R0);