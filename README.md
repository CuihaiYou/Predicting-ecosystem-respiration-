## Predicting-ecosystem-respiration
    This repository contains R scripts developed by Cuihai You for analyzing ecosystem respiration (ER) data in support of a publication in the Journal of Advances in Modeling Earth Systems (JAMES). The code performs four main analytical tasks: MCMC chain convergence diagnostics, slope comparison, ANOVA with post-hoc tests, and principal component analysis (PCA).
     Script Descriptions
   ## linear model
      main_linear.m – The main program for model execution
      Generate_C.m – Script for generating model parameters.
      soil_scalar_TECO.m – Calculates temperature and moisture scalars.
      solve_forward.m – Subroutine for solving the forward model
   ## microbial model
      main_microbial_model.m – The main program for model execution
      Generate_C.m – Script for generating model parameters.
      soil_scalar_TECO.m – Calculates temperature and moisture scalars.
      micro_process_new.m – Subroutine for solving the forward model
     
   ## Gelman-Rubin Convergence Diagnostic 
      This script assesses MCMC chain convergence for parameter estimation
      Reads four parallel MCMC chains from an Excel file (para_chain.xlsx)
      Converts each chain to an mcmc object using the coda package
      Combines chains into an mcmc.list object
      Calculates the Gelman-Rubin diagnostic (R-hat) to evaluate convergence
      Convergence is indicated by R-hat values close to 1 (typically <1.05)
      Usage Notes: 
      Replace **** with actual file path
      Requires the coda package for MCMC diagnostics

   ## Slope Comparison Analysis 
      This script evaluates model performance by comparing simulated vs. measured ER values:
      Reads data from Excel (your_data.xlsx)
      Fits a linear regression model (forced through origin)
      Generates model summary statistics
      Calculates 95% confidence intervals for the slope
      A slope not significantly different from 1 indicates good model performance
      Key Outputs:
      Regression coefficients and p-values
      Confidence interval for the slope parameter

   ## ANOVA and Multiple Comparisons 
      This script analyzes parameter differences across hydrological regimes:
      Reads parameter data (para_normal_dry_wet_linear.xlsx)
      Performs one-way ANOVA to test for treatment effects
      Conducts Tukey's HSD post-hoc tests for pairwise comparisons
      Identifies significant differences between normal, dry, and wet conditions
      Statistical Outputs:
      ANOVA table (F-statistic and p-value)
      Tukey test results showing pairwise differences

  ## Principal Component Analysis 
       This script explores parameter patterns across hydrological regimes:
       Reads parameter data (para_normal_dry_wet_linear.xlsx)
       Performs PCA using FactoMineR package
       Generates biplot visualization with fviz_pca_biplot
       Includes treatment grouping (normal/dry/wet) and 95% confidence ellipses
       Saves high-quality JPEG output (PCA_linear.jpeg)
       Visualization Features:
       Color-coded by treatment (red=normal, black=dry, blue=wet)
       Customized theme with publication-quality formatting
       Vector arrows showing parameter contributions
   ## Requirements
       R (version 3.6+ recommended)
       Required packages:
       readxl
       coda
       FactoMineR
       factoextra
       ggplot2


