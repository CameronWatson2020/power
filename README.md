# Power Analysis for PGS–Cognition Studies

This script performs power calculations and effect size estimations for polygenic score (PGS) associations with cognition across multiple datasets.

### Includes:
- Minimum detectable R² for UK Biobank, GAP/EUGEI, and MEG studies
- Sample size vs R² plots for different numbers of predictors
- Power calculations for connectivity models and SEM-based mediation

### Requirements:
```r
library(pwr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(lavaan)
library(semPower)

