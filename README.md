Code and Data for "Temperature and resources interact to affect transmission via host foraging rate and susceptibility" by Daniel C. Suh, Katie Schroeder, and Alexander T. Strauss

# Overview:

 - BUILD: instructions for rerunning code to generate data and figures
 - SETUP: descriptions of setup scripts
 - DATA: descriptions of raw and processed data
 - SCRIPTS: descriptions of scripts
 - Full Directory Structure: file tree
 - Session Info: R version, OS, and package versions used

# BUILD:


To reproduce analysis and figures, open temp_rsc_inter.Rproj and run _build.R. You must work within the .Rproj and have the 'here' package installed in order for all filepaths to work appropriately.
_build.R runs all scripts in sequence to generate all data and figures.

Scripts 00 through 05 perform all data cleaning and analysis.
Scripts 06 through 09 generate figures and statistical results.


# SETUP:


Necessary packages and some functions are specified in scripts in the "base" folder

├── base  
│   ├── fns.R     #custom functions  
│   ├── setup.R   #packages and plot theme  
│   └── src.R     #source fns.R and setup.R  


# FIGURES:


Figure components and tables are in the "figures" folder. Manuscript versions of figures with powerpoint file used for formatting are in the subfolder "manuscript_versions".

├── figures  
│   ├── 01_f_data.png   
│   ├── 01_length_data.png  
│   ├── 01_prev_data.png  
│   ├── 02_algae_plots.png  
│   ├── 02_clearance_plots.png  
│   ├── 03_inf_spores_plots.png  
│   ├── 03_prev_plots.png  
│   ├── 03_susc_plots.png  
│   ├── 04_handling.png  
│   ├── 04_param_est.png  
│   ├── 05_fora_ci_rsc.png  
│   ├── 05_fora_ci_temp.png  
│   ├── 05_susc_ci_rsc.png  
│   ├── 05_susc_ci_temp.png  
│   ├── boot_conf_int.docx            #bootstrapped CI for model 2E  
│   ├── combined_table_results.docx   #AIC table for infection model competition  
│   ├── fora_table_results.docx       #AIC table for foraging model competition  
│   ├── manuscript_versions  
│   │   ├── 01_data.png               #figure 1  
│   │   ├── 02_foraging.png           #figure 2  
│   │   ├── 03_infection.png          #figure 3  
│   │   ├── 04_handling_params.png    #figure 4  
│   │   ├── 05_bootstrap_f_u.png      #figure 5  
│   │   └── final_figures.pptx        #ppt for figure formatting  
│   ├── supp_body_length.png          #supplemental figure S1  
│   ├── supp_ci_est_combine.png       #supplemental figure S3  
│   └── supp_model_results.png        #supplemental figure S2  


# DATA:


Raw data files include:  
├── raw_data  
│   ├── day5_length.csv #length measurements used for infection models  
│   ├── foraging.csv    #foraging assay data  
│   ├── infection.csv   #infection assay data  
│   └── metadata.xlsx   #metadata for raw data files  

  day5_length.csv
  A representative subset of individuals were sacrificed to obtain length measurements on day 5 of the infection assay
  Needed for estimating body lengths of hosts at time of exposure
  Length measurements must be converted to mm depending on magnification level
  
  foraging.csv
  Raw data from foraging assay
  Needed for body length of hosts, duration of trials, and final fluorescence values
  Length measurements must be converted to mm depending on magnification level
  
  infection.csv
  Raw data from infection assay
  Needed for treatment conditions and infection status of all exposed hosts


Processed data files include:  
├── processed_data  
│   ├── foraging.rds                    #summary stats for foraging assay data  
│   ├── foraging_raw.rds                #cleaned foraging assay data  
│   ├── length.rds                      #summary stats for final lengths from infection assay  
│   ├── life_length_summ.rds            #summary stats for day 5 lengths from infection assay  
│   ├── m2E_bootstrap_quantiles.rds     #quantiles for bootstraps  
│   ├── m2E_bootstraps.rds              #raw bootstraps for best infection model  
│   ├── mle                             #maximum likelihood estimates
│   │   ├── m1_combined_fit.rds         #combined_fit are for infection models  
│   │   ├── m1_f_fit.rds                #f_fit are for foraging models  
│   │   ├── m2_combined_fit.rds  
│   │   ├── m2_f_fit.rds  
│   │   ├── m3_combined_fit.rds  
│   │   ├── m3_f_fit.rds  
│   │   ├── m4_combined_fit.rds  
│   │   ├── m4_f_fit.rds  
│   │   ├── m5_body_size_fit.rds        #size-dependent infection model (supp)  
│   │   ├── m5_combined_fit.rds  
│   │   └── m5_f_fit.rds  
│   ├── prevalence.rds                  #summary statistics for prevalence data  
│   ├── seq_data                        #simulated data w/ CI estimates  
│   │   ├── bootstrap_results.rds  
│   │   ├── f_ci_bootstrap.rds          
│   │   ├── foraging_rate_fit_data.rds  
│   │   ├── h_ci_bootstrap.rds          
│   │   ├── infection_fit_data.rds      
│   │   ├── supp_body_size_fit_data.rds   
│   │   └── u_ci_bootstrap.rds          
│   ├── supp_bootstrap_quantiles.rds    #quantiles for supp model bootstraps  
│   └── supp_bootstraps.rds             #raw bootstraps for supplemental model  


# SCRIPTS:


See below for a brief overview of each script:  

├── 00_data_processing.R  
├── 01_foraging.R  
├── 01a_foraging_mle.R  
├── 02_foraging_data.R  
├── 03_infection.R  
├── 03a_infection_mle.R  
├── 04_infection_data.R  
├── 05_bootstrap.R  
├── 06_tables_foraging.R  
├── 07_tables_infection.R  
├── 08_figures.R  
├── 09_stats.R  


  00_data_processing.R  
  Converts raw data into usable format for future analyses  
  Outputted data are saved in the "processed_data" folder  
  
  01_foraging.R  
  Defines model functions for foraging rate analysis  
  
  01a_foraging_mle.R  
  Estimates model parameters via maximum likelihood  
  Point estimates are saved in the "mle" folder under "processed_data"  
  This scripts sources 01_foraging.R  
  
  02_foraging_data.R  
  Simulates data using ML estimates   
  Data ares saved in the "seq_data" folder under "processed_data"  
  This scripts sources 01_foraging.R  
  
  03_infection.R  
  Defines model functions for foraging and transmission  
  
  03a_foraging.R  
  Estimates model parameters via maximum likelihood  
  Point estimates are saved in the "mle" folder under "processed_data"  
  This script sources 01_foraging.R  
  
  04_infection_data.R  
  Simulates data using ML estimates   
  Data ares saved in the "seq_data" folder under "processed_data"  
  This scripts sources 03_infection.R  
  
  05_bootstrap.R  
  Defines function and generates bootstrapped confidence interval estimates  
  Simulates confidence intervals for handling time, foraging rate, and per parasite susceptibility  
  Generates Table S3 in supplementary materials  
  This scripts sources 03_infection.R  
  
  06_tables_foraging.R  
  Generates AIC table for foraging models  
  Table 1A in manuscript  
  
  07_tables_infection.R  
  Generates AIC table for transmission models  
  Table 1B in manuscript  
  This scripts sources 03_infection.R  
  
  08_figures.R  
  Generates figure components for constructing manuscript figures  
  Final manuscript versions of figures were formatted using Microsoft Powerpoint  
  Figures 1-5 in manuscript  
  
  09_stats.R  
  Performs traditional statistics  
  Reported in text and in figure 1 in manuscript  


Supplemental Analysis and Figures  

├── S0_body_length_fig.R  
├── S1_body_size_fit.R  
├── S2_seq_data.R  
├── S3_bootstrap.R  

  S0_body_length_fig.R  
  Compares body lengths between foraging and infection assays  
  Generates supplementary figure S1  

  S1_body_size_fit.R  
  Estimates model parameters via maximum likelihood for additional model (2F)  
  This scripts sources 03_infection.R  

  S2_seq_data.R  
  Simulates data using ML estimates  
  Data ares saved in the "seq_data" folder under "processed_data"  
  This scripts sources 03_infection.R  
  Generates supplementary figure S2  

  S3_bootstrap.R  
  Generates bootstrapped confidence intervals for model 2F  
  This scripts sources 03_infection.R  
  Generates supplementary figure S3  




# Full Directory Structure:


├── 00_data_processing.R  
├── 01_foraging.R  
├── 01a_foraging_mle.R  
├── 02_foraging_data.R  
├── 03_infection.R  
├── 03a_infection_mle.R  
├── 04_infection_data.R  
├── 05_bootstrap.R  
├── 06_tables_foraging.R  
├── 07_tables_infection.R  
├── 08_figures.R  
├── 09_stats.R  
├── README.md  
├── S0_body_length_fig.R  
├── S1_body_size_fit.R  
├── S2_seq_data.R  
├── S3_bootstrap.R  
├── _build.R  
├── base  
│   ├── fns.R  
│   ├── setup.R  
│   └── src.R  
├── figures  
│   ├── 01_f_data.png  
│   ├── 01_length_data.png  
│   ├── 01_prev_data.png    
│   ├── 02_algae_plots.png  
│   ├── 02_clearance_plots.png  
│   ├── 03_inf_spores_plots.png  
│   ├── 03_prev_plots.png  
│   ├── 03_susc_plots.png  
│   ├── 04_handling.png  
│   ├── 04_param_est.png  
│   ├── 05_fora_ci_rsc.png  
│   ├── 05_fora_ci_temp.png  
│   ├── 05_susc_ci_rsc.png  
│   ├── 05_susc_ci_temp.png  
│   ├── boot_conf_int.docx  
│   ├── combined_table_results.docx  
│   ├── fora_table_results.docx  
│   ├── manuscript_versions  
│   │   ├── 01_data.png  
│   │   ├── 02_foraging.png  
│   │   ├── 03_infection.png  
│   │   ├── 04_handling_params.png  
│   │   ├── 05_bootstrap_f_u.png  
│   │   └── final_figures.pptx  
│   ├── supp_body_length.png  
│   ├── supp_ci_est_combine.png  
│   └── supp_model_results.png  
├── processed_data  
│   ├── foraging.rds  
│   ├── foraging_raw.rds  
│   ├── length.rds  
│   ├── life_length_summ.rds  
│   ├── m2E_bootstrap_quantiles.rds  
│   ├── m2E_bootstraps.rds  
│   ├── mle  
│   │   ├── m1_combined_fit.rds  
│   │   ├── m1_f_fit.rds  
│   │   ├── m2_combined_fit.rds  
│   │   ├── m2_f_fit.rds  
│   │   ├── m3_combined_fit.rds  
│   │   ├── m3_f_fit.rds  
│   │   ├── m4_combined_fit.rds  
│   │   ├── m4_f_fit.rds  
│   │   ├── m5_body_size_fit.rds  
│   │   ├── m5_combined_fit.rds  
│   │   └── m5_f_fit.rds  
│   ├── prevalence.rds  
│   ├── seq_data  
│   │   ├── bootstrap_results.rds  
│   │   ├── f_ci_bootstrap.rds  
│   │   ├── foraging_rate_fit_data.rds  
│   │   ├── h_ci_bootstrap.rds  
│   │   ├── infection_fit_data.rds  
│   │   ├── supp_body_size_fit_data.rds  
│   │   └── u_ci_bootstrap.rds  
│   ├── supp_bootstrap_quantiles.rds  
│   └── supp_bootstraps.rds  
├── raw_data  
│   ├── day5_length.csv  
│   ├── foraging.csv  
│   ├── infection.csv  
│   └── metadata.xlsx  
└── temp_rsc_inter.Rproj  


# Session Info:  


R version 4.3.1 (2023-06-16)  
Platform: aarch64-apple-darwin20 (64-bit)  
Running under: macOS Sonoma 14.1.1  

Matrix products: default  
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib   
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0  

locale:  
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8  

time zone: America/New_York  
tzcode source: internal  
  
attached base packages:  
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2          readr_2.1.4          tidyr_1.3.1          tibble_3.2.1         tidyverse_2.0.0     
 [9] patchwork_1.3.0.9000 egg_0.4.5            gridExtra_2.3        AICcmodavg_2.3-3     flextable_0.9.6      ggnewscale_0.4.9     bbmle_1.0.25         deSolve_1.36        
[17] lhs_1.1.6            rstatix_0.7.2        ggpubr_0.6.0         ggplot2_3.5.1        epitools_0.5-10.1    lubridate_1.9.4      magrittr_2.0.3       here_1.0.1          

loaded via a namespace (and not attached):
 [1] pbapply_1.7-2           rlang_1.1.4             compiler_4.3.1          systemfonts_1.0.4       vctrs_0.6.5             httpcode_0.3.0          pkgconfig_2.0.3        
 [8] crayon_1.5.3            fastmap_1.1.1           backports_1.5.0         ellipsis_0.3.2          labeling_0.4.3          promises_1.2.1          rmarkdown_2.29         
[15] tzdb_0.4.0              ragg_1.2.5              bit_4.0.5               xfun_0.49               jsonlite_1.8.7          later_1.3.1             uuid_1.1-0             
[22] VGAM_1.1-9              broom_1.0.5             parallel_4.3.1          R6_2.5.1                stringi_1.8.4           car_3.1-2               numDeriv_2016.8-1.1    
[29] Rcpp_1.0.11             knitr_1.43              httpuv_1.6.11           Matrix_1.5-4.1          splines_4.3.1           timechange_0.3.0        tidyselect_1.2.1       
[36] rstudioapi_0.15.0       abind_1.4-8             curl_6.1.0              lattice_0.21-8          shiny_1.7.5             withr_3.0.2             askpass_1.1            
[43] evaluate_0.21           survival_3.5-5          zip_2.3.1               xml2_1.3.5              pillar_1.10.1           carData_3.0-5           generics_0.1.3         
[50] vroom_1.6.3             rprojroot_2.0.3         hms_1.1.3               munsell_0.5.1           scales_1.3.0            xtable_1.8-4            glue_1.8.0             
[57] unmarked_1.3.2          gdtools_0.3.7           tools_4.3.1             gfonts_0.2.0            data.table_1.14.8       ggsignif_0.6.4          mvtnorm_1.3-2          
[64] grid_4.3.1              bdsmatrix_1.3-6         colorspace_2.1-1        nlme_3.1-162            cli_3.6.3               textshaping_0.3.6       officer_0.6.6          
[71] fontBitstreamVera_0.1.1 gtable_0.3.6            digest_0.6.33           fontquiver_0.2.1        crul_1.4.2              farver_2.1.2            htmltools_0.5.6        
[78] lifecycle_1.0.4         mime_0.12               bit64_4.0.5             fontLiberation_0.1.0    openssl_2.1.0           MASS_7.3-60  

