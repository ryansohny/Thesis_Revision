library(survival)
library(survminer)

"""
> sessionInfo()
R version 4.1.3 (2022-03-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS release 6.10 (Final)

Matrix products: default
BLAS/LAPACK: /mnt/data/Projects/phenomata/99.Tools/Anaconda3/envs/R/lib/libopenblasp-r0.3.20.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] survminer_0.4.9 ggpubr_0.4.0    ggplot2_3.3.6   survival_3.3-1 

loaded via a namespace (and not attached):
 [1] pillar_1.7.0      compiler_4.1.3    tools_4.1.3       digest_0.6.29    
 [5] lifecycle_1.0.1   tibble_3.1.7      gtable_0.3.0      lattice_0.20-45  
 [9] pkgconfig_2.0.3   rlang_1.0.2       Matrix_1.4-1      cli_3.3.0        
[13] xfun_0.31         gridExtra_2.3     knitr_1.39        withr_2.5.0      
[17] dplyr_1.0.9       generics_0.1.2    vctrs_0.4.1       survMisc_0.5.6   
[21] grid_4.1.3        tidyselect_1.1.2  data.table_1.14.2 glue_1.6.2       
[25] KMsurv_0.1-5      R6_2.5.1          rstatix_0.7.0     km.ci_0.5-6      
[29] fansi_1.0.3       carData_3.0-5     farver_2.1.0      purrr_0.3.4      
[33] tidyr_1.2.0       car_3.1-0         magrittr_2.0.3    scales_1.2.0     
[37] backports_1.4.1   ellipsis_0.3.2    splines_4.1.3     abind_1.4-5      
[41] xtable_1.8-4      colorspace_2.0-3  ggsignif_0.6.3    labeling_0.4.2   
[45] utf8_1.2.2        munsell_0.5.0     broom_0.8.0       crayon_1.5.1     
[49] zoo_1.8-10
"""

tcgaov <- read.csv("TCGA_clinical_2022.csv")
# Column Info
# status ==> '2': DECEASED, '1': LIVING
# emt_s ==> '1': EMT-high, '2': EMT-low
# emt_g ==> '1': Mesenchymal, '2': 'Epithelial'
# emt_c ==> '1': EMT(+), '2': EMT(-)

# For EMT index-based classification (A Method from Sohn et al.)
fit <- survfit(Surv(time, status) ~ emt_s, data = tcgaov)
ggsurvplot(fit, conf.int = TRUE, pval = TRUE, palette = c('#990000', '#000066'))

# For EMT score-based classification (A Method from Guo et al.)
fit <- survfit(Surv(time, status) ~ emt_g, data = tcgaov)
ggsurvplot(fit, conf.int = TRUE, pval = TRUE, palette = c('#990000', '#000066'))

# For EMT index-based classification (A Method from Cristescu et al.)
fit <- survfit(Surv(time, status) ~ emt_c, data = tcgaov)
ggsurvplot(fit, conf.int = TRUE, pval = TRUE, palette = c('#990000', '#000066'))
