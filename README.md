# MMB-117_Environmental_microbiology
**Master's Programme in Microbiology and Biotechnology — Advanced Level**
**5 credits | University of Helsinki**

Course was last time held in January 2025. This is the data analysis part.

To run the analysis script MMB-117.R, you will need to download that and the files seqtab.rds and taxa.rds, open the file MMB-117.R in your R and follow the instructions in the scrip.

The script is for analyzing soil microbiomes in four sites and what environmental variables explain the microbiomes 

# Soil Microbiome Analysis Pipeline
---

## About the Course

MMB-117 is a research-oriented field and laboratory course where students experience the full arc of a microbiome study — from hypothesis to analysis. The course unfolds in four phases:

1. **Research planning** — Students attend a theory lecture on scientific thinking and hypothesis formulation, then participate in a facilitated workshop to define research questions and design a sampling plan. Link to theory lecture slides: https://docs.google.com/presentation/d/15MKmCsdzrR8Ai7pkKF3dFHeQsE79-Jq2/edit?slide=id.p1#slide=id.p1 
2. **Field sampling** — The following day, the class collects soil samples together according to the student-designed plan.
3. **Laboratory work (~2 weeks)** — Students perform standard soil analyses (gravimetric water content, pH, soil organic matter via loss-on-ignition) and extract DNA from each sample. The DNA extracts are sent for **16S rRNA gene amplicon sequencing**.
4. **Bioinformatic and statistical analysis** — This repository contains the R pipeline used to analyse the sequencing data, using the laboratory measurements and sampling sites as metadata.

---

## Repository Contents

```
├── README.md
├── MMB117_analysis.R          # Main analysis script (documented below)
├── MMB117_Sample_sheet_2025.txt   # Sample metadata 
├── seqtab.rds                 # ASV count table from DADA2 
└── taxa.rds                   # Taxonomy table from DADA2

```
---

## Analysis Overview

### 1. Data Import and Preparation

- Load ASV count table (`seqtab.rds`) and taxonomy table (`taxa.rds`) produced by DADA2
- Rename ASVs sequentially (`ASV1`, `ASV2`, …) for readability
- Apply **Total Sum Scaling (TSS)** to convert counts to relative abundances
- Remove the negative sequencing control from the ASV table
- Align sample order between metadata and ASV table
- Calculate **Shannon's alpha diversity index** per sample and add to metadata

### 2. Data Exploration

Following the protocol of [Zuur et al. (2010)](https://doi.org/10.1111/j.2041-210X.2009.00001.x):

| Exploration step | Variables checked |
|---|---|
| Outliers in X (metadata) | WetWeight, DryWeight, DNAWeight, Moisture, pH (water & CaCl₂), SOM |
| Outliers in Y (microbiome) | Shannon diversity, NMDS ordination |
| Homogeneity of variances | Boxplots of all metadata variables by site |
| Normality | Discussed theoretically; residuals to be checked after modelling |
| Zero inflation | Visualised in the relative abundance matrix |
| Collinearity in X | Pairwise scatter plots of soil variables |
| Relationships Y ~ X | Scatter plots of diversity against soil variables |
| Interactions | Visualised with faceted plots (sample size too small to model) |

**Key findings from exploration:**
- pH (water and CaCl₂) are collinear — only pH_CaCl₂ is carried forward
- Moisture shows high heteroscedasticity across sites (likely influenced by snow at sampling time)
- No extreme outliers warranting removal; negative control behaves as expected

### 3. Ordination (NMDS)

- Non-metric multidimensional scaling (NMDS) on TSS-normalised Bray-Curtis distances (`vegan::metaMDS`)
- **Environmental fitting** of pH_CaCl₂ and Shannon diversity as linear vectors (`vegan::envfit`)
  - pH_CaCl₂ explains **~92%** of community variance (p < 0.001)
  - Shannon diversity explains **~82%** of community variance (p < 0.001)
- **Surface fitting** (GAM) of soil organic matter (SOM) and gravimetric water content (GWC) via `vegan::ordisurf`
  - SOM: R² = 0.79, deviance explained = 86% (p < 0.001)
  - GWC: R² = 0.67, deviance explained = 76% (p < 0.001)

### 4. PERMANOVA

Community composition tested with `vegan::adonis2` (9999 permutations, sequential terms):

| Model | Key result |
|---|---|
| Site + pH + SOM + Moisture | Only **Site** significant (R² = 0.58, p < 0.001) |
| pH + SOM + Moisture (no Site) | **pH** (R² = 0.31, p < 0.001) and **Moisture** (p = 0.009) significant |
| pH alone | R² = 0.27, p < 0.001 |
| Site alone | R² = 0.57, p < 0.001 |

Site and pH are collinear (sites differ in pH), so they cannot be interpreted independently.

### 5. Alpha Diversity

- Shannon diversity compared across sites using **t-tests** with `ggpubr::stat_compare_means`, Forest as reference group
- One-way ANOVA + Tukey HSD as a complementary test
- **Gas station** and **Park** had significantly higher diversity than Forest (p < 0.001)

---

## R Packages Required

```r
install.packages(c(
  "reshape2", "dplyr", "tidyr", "MASS", "vegan",
  "ggplot2", "ggthemes", "colorspace", "colorRamps",
  "RColorBrewer", "gridExtra", "fitdistrplus", "logspline",
  "car", "lme4", "scales", "viridis", "ggrepel",
  "pheatmap", "ggpubr"
))
```

---

## Discussion Points for Students

- Why is **Site** the strongest predictor — does it represent land-use history, pH, or both?
- Shannon diversity is **highest at the gas station** — what ecological processes could explain this?
- Why can't we include both Site and pH in the same PERMANOVA model?
- What are the limitations of TSS normalisation, and when would CLR transformation be preferable?
- How does the amount of snow at sampling time complicate the interpretation of gravimetric water content?

---

## References

- Zuur, A.F., Ieno, E.N. and Elphick, C.S. (2010). A protocol for data exploration to avoid common statistical problems. *Methods in Ecology and Evolution*, 1: 3–14. https://doi.org/10.1111/j.2041-210X.2009.00001.x
- Oksanen J et al. (2024). *vegan: Community Ecology Package*. R package. https://CRAN.R-project.org/package=vegan
- Callahan BJ et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13: 581–583.

---

## Contact

For questions about the course or pipeline, please email: johanna.muurinen@helsinki.fi
