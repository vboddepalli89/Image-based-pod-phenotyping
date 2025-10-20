**AI-assisted Image-Based Phenotyping Reveals Genetic Architecture of Pod Traits in Mungbean**

#Overview
This repository accompanies the manuscript “AI-assisted Image-Based Phenotyping Reveals Genetic Architecture of Pod Traits in Mungbean.”
It provides all R scripts and supplementary data used for image-based trait extraction, statistical analyses, and genomic association mapping.
The project integrates high-throughput phenotyping and genome-wide analyses to dissect the genetic basis of pod morphology and yield-related traits in mungbean (Vigna radiata).
| File                               | Description                                                                                             |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------- |
| **01_descriptive_statistics.R**    | Performs descriptive statistical analyses and trait distribution visualization.                         |
| **02_methods_comparison_fig2.Rmd** | Compares manual and image-based trait measurements (Figure 2).                                          |
| **03_trait_blups.Rmd**             | Computes Best Linear Unbiased Predictors (BLUPs) using mixed linear models.                             |
| **04_gwas_gapit_pipeline.Rmd**     | Implements Genome-Wide Association Study (GWAS) using GAPIT with multiple models (BLINK, FarmCPU, MLM). |
| **05_gwas_sven_method.R**          | Executes the SVEN GWAS method for comparison with GAPIT results.                                        |
| **06_rrblup_genomic_prediction.R** | Conducts genomic prediction (GP) using rrBLUP and evaluates predictive accuracy.                        |
| **s02_venn_diagram_figS2.R**       | Generates Venn diagrams for overlapping significant SNPs (Supplementary Figure S2).                     |

Methodological Summary
Pod and seed traits were extracted through an AI-assisted image analysis pipeline developed for automated phenotyping.
Phenotypic data were analyzed using mixed linear models to estimate BLUPs and heritability.
Genome-wide association studies were performed using GAPIT and SVEN to identify loci associated with pod traits.
Genomic prediction analyses were conducted with rrBLUP to assess the accuracy of genomic selection.
All analyses were performed in the R environment (v4.2.0 or later).

All analyses were performed using R (version 4.2.0 or later).
Key packages and versions used include:
  GAPIT3 (version 3.1.0) – genome-wide association analysis
  rrBLUP (version 4.6.1) – genomic prediction
  sommer (version 4.3.3) – mixed model analysis
  lme4 (version 1.1-34) – linear mixed-effects models
  data.table (version 1.14.8) – high-performance data manipulation
  ggplot2 (version 3.4.2) – visualization
  dplyr (version 1.1.2) – data wrangling
