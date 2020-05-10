%let coviddir=/folders/myfolders/COVID19_SAS;
%include "&coviddir/src/calchpd.sas";
libname covid "&coviddir/saslib";
options compress=char;