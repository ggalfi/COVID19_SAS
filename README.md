# COVID19_SAS

My intentions was two-folded when I went through this analysis.

First: I wanted to get answer for a frequently arising question: how many people have COVID-19 in a certain country/region. Having this number is a must for a lot of calculatons, for example to estimate somehow the mortality rate. But the problem is that the officially available data is highly unreliable becasue of a few good reasons. So I tried to find the answer with Bayesian analysis. What even surprised me, that pretty stable and consistent results were came out from this analysis.
See the [report](https://github.com/ggalfi/COVID19_SAS/raw/master/doc/covid19mcmc.pdf) I have written on this topic.

Second: This analysis is a very good example of usefullnes of the Markov Chain Monte Carlo method which I've relied on during the analysis. So [grab somewhere a SAS University Edition](https://www.sas.com/en_us/software/university-edition/download-software.html) (don't worry, it's free) download this repository which contains all the data and sas code I've used in this report and you are ready to jump into the interesting world of Bayesian analysis.

To make SAS codes work, clone this repository into the shared folder of SAS UE. Then open initlib.sas in SASStudio and check that "coviddir" macro variable points to this repository. Change it when needed. Run the whole SAS code and then all set for further analysis. You may run exploration.sas to see the time series which we want to apply MCMC on (it is optional). By running mcmc.sas the base datasets which my report relies on are created, e.g. a posterior sample of parameters. Feel free to play with the settings of MCMC (prior distributions, sample size, optimization method for proposal function, etc.). After having the MCMC result you can run post_mcmc_analysis.sas which does all the calculations I used in my report.


