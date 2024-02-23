# Data and R scripts for the manuscript 'Modelling Diffusion of Innovation Curves using Radiocarbon Data'

This repository contains data and scripts used in the manuscript:

Crema, E.R., Bloxam, A., Stevens, C.J., Vander Linden, M. (2023). Modelling Diffusion of Innovation Curves using Radiocarbon Data

The repository is organised into the following six main directories:

- **analysis** ... Contains R scripts for analysing the empirical data sets.
- **data** ... Contains raw data and pre-processing scripts for the empirical case studies.
- **figures_and_tables** ... Contains R scripts for generating figures and tables included in the manuscript and the supplementary data.
- **results** ... Contains R images files of outputs the Bayesian Analysis.
- **sim** ... Contains R scripts for the creation and the analyses of the simulated datasets.
- **src** ... Contains custom R functions.

## Dataset
Three empirical case studies case studies are examined:

_Case Study 1a_ :
 - Diffusion of Rice agriculture in Prehistoric Japan 4000-1700 cal BP (`jpdata.RData`).
   
_Case Study 1b_ :
 - Diffusion of Farming in Prehistoric Britain, 7000-3000 cal BP (`gbdata.RData`).
   
_Case Study 2_:
 - Diffusion cycles of Cremation and Inhumation in Prehistoric Britain, 5500-2200 cal BP `burialdata.RData`. 

Empirical datasets for case studies I and II were obtained from the following resources and publications:
- Crema, E. R., Stevens, C. J., & Shoda, S. (2022). Bayesian analyses of direct radiocarbon dates reveal geographic variations in the rate of rice farming dispersal in prehistoric Japan. Science Advances, 8(38), eadc9171. (File: `/Data/R14CDB.csv`, retrieved from https://github.com/ercrema/yayoi_rice_dispersal)
- [Database of radiocarbon dates published in Japanese Archaeological Site Reports](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd_en/db_param ) (File: `data/c14db_1.1.0.csv`)
- Bevan, A., Colledge, S., Fuller, D., Fyfe, R., Shennan, S. Stevens, C.(2017). Holocene Fluctuations in Human Population Demonstrate Repeated Links to Food Production and Climate. Proceedings of the National Academy of Sciences 114(49) : E10524–31. https://doi.org/10.1073/pnas.1709190114. (retrieved directly [from UCL Discovery Repository](https://discovery.ucl.ac.uk/id/eprint/10025178/4/Bevan_gbie14Csub.zip) in the script `/data/gb_clean.R`)
  
## Simulation
Simulated datasets with sigmoid growth (`simdata1a.RData` and `simdata1b.RData`) were generated using the script `simulate_sigmoid.R` using sample sizes and number of locations equal to the two case studies (`jpdata.RData` and `gbdata.RData`). The script `simulate_icar.R` generates instead an arbitrary non-linear diffusion dynamic (`simdata2.RData`) using the same sample size as the burial dataset (`burialdata.RData`).   

## Bayesian Analysis
R scripts for core Bayesian analyses on the empirical case studies and simulated datasets are stored respectively in the `analysis` and `sim` directories. Model fitting is executed via Nimble probabilistic programming language and requires around 24 to be completed. In the case of the simulated datasets and case study 2 the output of each script (`fit_sim1a.R`,`fit_sim1b.R`, `fit_sim2.R` and `burial_icar.R`) are stored within individual R image files (`post_sim1a.RData`, `post_sim1b.RData`,  `post_icar_sim2.RData`  and `post_icar_burial.RData` ). In the case of the sigmoid models of case studies 1a and 1b, the core scripts `japan_abot.R` and `britain_abot.R` produces two R image files respectively, one containing the posterior of the parameters (`post_jp_abot.RData` and `post_gb_abot.RData`) and the other containing 1,000 sets of simulated data points from random samples of the posterior parameters (`ppcheck_jp_abot.RData` and `ppcheck_gb_abot.RData`). The latter is then read by a dedicated R script (`post_check_jp_abot.R` and `post_check_gb_abot.R`) which executes the posterior predictive checks (i.e. by generating envelops of simulated proportion summed probability distribution of radiocarbon dates) and stores the output into another R image file (`ppc_jp_abot.RData` and `ppc_gb_abot.RData`).  

## File Structure

### analysis
 - `britain_abot.R` ... fits a hierarchical sigmoid model on `gbdata.RData`; outputs `post_gb_abot.RData` and `ppc_gb_abot.RData`.    
 - `burial_icar.R` ... fits an ICAR model on `burialdata.RData`; outputs `post_icar_burial.RData`
 - `japan_abot.R` ... fits a hierarchical sigmoid model on `jpdata.RData`; outputs `post_jp_abot.RData` and `ppc_jp_abot.RData`. 
 - `post_check_gb_abot.R` ... runs posterior predictive checks on `ppc_gb_abot.RData`; outputs `ppcheck_gb_abot.RData`
 - `post_check_jp_abot.R` ... runs posterior predictive checks on `ppc_jp_abot.RData`; outputs `ppcheck_jp_abot.RData`
### data
 - `burial_clean.R` ... reads and cleans `burialdates.csv`; outputs `burialdata.RData`  
 - `burialdata.RData` ... burial data for case study 2
 - `gb_clean.R` ... downloads and processes data for case study I, Britain; outputs `gbdata.RData` 
 - `gbdata.RData` ... archaeobotanical data for case study 1a, Britain.
 - `jp_clean.R` ... reads, cleans, and combines `R14CDB.csv`, `c14db_1.1.0.csv`, and `Taxa_Edible_Classifications.csv` ; outputs `jpdata.RData` 
 - `jpdata.RData` ... archaeobotanical data for case study 1b, Japan.
 -  _raw_
     - `R14CDB.csv` ... rice radiocarbon database, obtained from https://github.com/ercrema/yayoi_rice_dispersal
     - `Taxa_Edible_Classifications.csv` ... lookup table to classify archaeobotanical samples into rice and wild nuts categories
     - `burialdates.csv` ... burial custom radiocarbon database
     - `c14db_1.1.0.csv` ... radiocarbon database of Japan, obtained from https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd_en/db_param

### figures_and_tables
  - `figure1.pdf` ~ `figure5.pdf` ... main figures in the manuscript.
  - `figureS1.pdf` ~ `figureS6.pdf` ... supplementary figures.
  - `figures_main.R` ... R scripts for generating main figures.
  - `figures_esm.R` ... R scripts for generating supplementary figures.
  - `table_main.R` ... R script for generating table 1
  - `table1.csv` ... table 1
### results 
  - `post_gb_abot.RData` ... posterior estimates for the hierarchical sigmoid model, case Study 1a, Britain. 
  - `post_icar_burial.RData`  ... posterior estimates for the ICAR model, case Study 2.
  - `post_jp_abot.RData`  ... posterior estimates for the hierarchical sigmoid model, case Study 1a, Britain.
  - `ppc_gb_abot.RData` ... posterior predictive estimates for the hierarchical sigmoid model, case Study 1a, Britain. 
  - `ppc_jp_abot.RData` ... posterior predictive estimates for the hierarchical sigmoid model, case Study 1b Japan. 
  - `ppcheck_gb_abot.RData` ... posterior predictive checks for the hierarchical sigmoid model, case Study 1a, Britain.  
  - `ppcheck_jp_abot.RData` ... posterior predictive checks for the hierarchical sigmoid model, case Study 1b, Japan.  
### sim
  - `fit_sim1a.R` ... fits a hierarchical sigmoid model on `simdata1a.RData`; outputs `post_sim1a.RData`.   
  - `fit_sim1b.R` ... fits a hierarchical sigmoid model on `simdata1b.RData`; outputs `post_sim1b.RData`.
  - `fit_sim2.R` ... fits an ICAR model on `simdata2.RData`; outputs `post_icar_sim2.RData`
  - `simulate_icar.R` ... simulates non-parametric cyclical diffusion data; outputs `simdata2.RData.
  - `simulate_sigmoid.R` ... simulates sigmoid diffusion data; outputs `simdata1a.RData` and `simdata1b.RData`
  -  _simdata_
     - `simdata1a.RData` ... simulated dataset 1a
     - `simdata1b.RData` ... simulated dataset 1b
     - `simdata2.RData` ... simulated dataset 2
  -  _simdata_
     - `post_sim1a.RData` ... posterior estimates for the hierarchical sigmoid model, simulated dataset 1a
     - `post_sim1b.RData` ... posterior estimates for the hierarchical sigmoid model, simulated dataset 1b
     - `post_icar_sim2.RData` ... posterior estimates for the ICAR model, simulated dataset 2
### src
  - `utility.R` ... Variety of utility functions including for basic calculations, posterior predictive checks, and plotting.

## R Session Info
```
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] truncnorm_1.0-9     here_1.0.1          coda_0.19-4        
 [4] latex2exp_0.9.6     RColorBrewer_1.1-3  dplyr_1.1.2        
 [7] nimbleCarbon_0.2.4  nimble_1.0.1        sf_1.0-14          
[10] rnaturalearth_0.3.3 rcarbon_1.5.0      

loaded via a namespace (and not attached):
 [1] xfun_0.40              spatstat.sparse_3.0-2  lattice_0.21-8        
 [4] numDeriv_2016.8-1.1    vctrs_0.6.2            tools_4.3.1           
 [7] doSNOW_1.0.20          spatstat.utils_3.0-3   generics_0.1.3        
[10] goftest_1.2-3          parallel_4.3.1         tibble_3.2.1          
[13] proxy_0.4-27           fansi_1.0.4            pkgconfig_2.0.3       
[16] spatstat_3.0-6         Matrix_1.6-0           KernSmooth_2.23-22    
[19] lifecycle_1.0.3        compiler_4.3.1         stringr_1.5.0         
[22] deldir_1.0-9           spatstat.linnet_3.1-1  codetools_0.2-19      
[25] spatstat.explore_3.2-1 snow_0.4-4             class_7.3-22          
[28] pracma_2.4.2           pillar_1.9.0           classInt_0.4-9        
[31] spatstat.model_3.2-4   iterators_1.0.14       rpart_4.1.19          
[34] abind_1.4-5            foreach_1.5.2          nlme_3.1-162          
[37] spatstat.geom_3.2-4    tidyselect_1.2.0       stringi_1.7.12        
[40] splines_4.3.1          rprojroot_2.0.3        polyclip_1.10-4       
[43] grid_4.3.1             cli_3.6.1              magrittr_2.0.3        
[46] utf8_1.2.3             e1071_1.7-13           spatstat.data_3.0-1   
[49] tensor_1.5             sp_2.0-0               httr_1.4.6            
[52] igraph_1.5.1           knitr_1.43             mgcv_1.9-0            
[55] rlang_1.1.1            Rcpp_1.0.11            spatstat.random_3.1-5 
[58] glue_1.6.2             DBI_1.1.3              jsonlite_1.8.5        
[61] R6_2.5.1               units_0.8-3  
```


## Funding
* ERC Starting Grant _Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan_ (ENCOUNTER) (Project N. 801953, PI: E. Crema).
* Philip Leverhulme Prize (#PLP-2019–304 Awarded to: E.Crema)
* NERC radiocarbon grant (#NF/2017/2/12, A. Bloxam, Awarded to: M.Parker Pearson)
* LACHP AHRC Studentship (Awarded to: A. Bloxam)
  
## Licence
CC-BY 3.0
