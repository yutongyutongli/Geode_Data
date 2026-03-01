# Geode Dataset: Preprocessed Data and Analysis Guide

## Overview

This repository contains a preprocessed dataset (`preprocessed_data.RData`) and comprehensive analysis code (`data_analysis_final.R`) for a study examining cognitive performance differences between typically developing (TD) and autism spectrum disorder (ASD) participants. The dataset includes behavioral data from a visual discrimination task involving counting flashes, along with demographic information and clinical assessments.

### Task Description

The data comes from a visual discrimination task where participants were presented with flashing stimuli on the left and right sides of the screen. Participants had to determine which side had more flashes and respond accordingly. The task included:

- **Stimulus presentation**: Flashes were presented in 10 temporal bins on each side
- **Response collection**: Participants made binary choices (left vs right)
- **Feedback**: Participants received feedback on their performance
- **Learning component**: The task included 200 trials to assess learning over time

The task was designed to assess visual processing, numerical discrimination, and learning abilities in different participant groups.

## Dataset Structure

### Main Data Object: `alldata`

The `preprocessed_data.RData` file contains a single data.table object called `alldata` with the following structure:

#### Core Variables
- **`username`**: Unique participant identifier
- **`subID`**: Subject ID (may include suffixes like "_2", "_3" for multiple sessions)
- **`game_version`**: Task version ("ft" = first time, "ft-2", "ft-3" for repeat sessions)
- **`source`**: Data source identifier
- **`group`**: Original group classification (includes "ASD" and "TD" variants)
- **`trial`**: Trial number within the task (1-200 for complete sessions)
- **`score`**: Binary score for each trial (0 or 1)
- **`rt`**: Response time in milliseconds
- **`perf`**: Overall performance (proportion correct across all trials)
- **`trainperf`**: Training performance
- **`nTrials`**: Number of completed trials

#### Task-Specific Variables
- **`choice`**: Participant's choice (typically 0 for left, 1 for right, or similar binary coding)
- **`dflash`**: Difference in number of flashes between left and right sides (can be negative or positive)
- **`tflash`**: Total number of flashes presented on both sides
- **`numflashleft`**: Number of flashes presented on the left side
- **`numflashright`**: Number of flashes presented on the right side
- **`winstay`**: Win-stay behavior indicator (1 if participant repeated previous choice after a win, 0 otherwise)
- **`loseswitch`**: Lose-switch behavior indicator (1 if participant switched choice after a loss, 0 otherwise)

#### Temporal Flash Variables (10 time bins)
- **`lbin1` through `lbin10`**: Binary indicators for flashes in each of 10 time bins on the left side
- **`rbin1` through `rbin10`**: Binary indicators for flashes in each of 10 time bins on the right side

#### Derived Temporal Variables
- **`l_early`**: Sum of flashes in early time bins on left side (lbin1 + lbin2 + lbin3)
- **`l_middle`**: Sum of flashes in middle time bins on left side (lbin4 + lbin5 + lbin6 + lbin7)
- **`l_late`**: Sum of flashes in late time bins on left side (lbin8 + lbin9 + lbin10)
- **`r_early`**: Sum of flashes in early time bins on right side (rbin1 + rbin2 + rbin3)
- **`r_middle`**: Sum of flashes in middle time bins on right side (rbin4 + rbin5 + rbin6 + rbin7)
- **`r_late`**: Sum of flashes in late time bins on right side (rbin8 + rbin9 + rbin10)

#### Demographic Variables
- **`gender`**: Participant gender ("M", "F", or other)
- **`age`**: Participant age in years
- **`race`**: Race classification (coded values that get mapped to categories)
- **`ethnicity`**: Ethnicity classification (coded values)

#### Clinical Assessment Variables
- **`vinabcstd`**: Vineland Adaptive Behavior Scales - 3rd Edition standard score
- **`srs2total`**: Social Responsiveness Scale - 2nd Edition total score
- **`bistotal`**: Behavior Rating Inventory of Executive Function total score
- **`aasp_low_reg_raw`**: Adolescent/Adult Sensory Profile - Low Registration raw score
- **`aasp_sen_seek_raw`**: Adolescent/Adult Sensory Profile - Sensory Seeking raw score
- **`aasp_sen_ses_raw`**: Adolescent/Adult Sensory Profile - Sensory Sensitivity raw score
- **`aasp_sen_avoid_raw`**: Adolescent/Adult Sensory Profile - Sensory Avoiding raw score

#### Derived Variables (Created in Analysis)
- **`binary_group`**: Simplified group classification ("ASD" vs "TD")
- **`group_verbal`**: Verbal ability classification based on VABS-3 scores
- **`srs_normal`**: Binary indicator for SRS-2 scores meeting inclusion criteria
- **`above_crit`**: Binary indicator for performance and trial completion criteria

## Data Analysis Code

### Prerequisites

Before running the analysis, ensure you have the following R packages installed:

```r
# Core data manipulation and visualization
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(tidyverse)

# Statistical analysis
library(Hmisc)
library(rstan)
library(zoo)
library(afex)
library(emmeans)
library(minpack.lm)
library(glmnet)
library(lme4)
library(broom.mixed)

# Utilities
library(latex2exp)
```

### Running the Analysis

1. **Set up your working directory**:
   ```r
   # Update this path to match your local directory structure
   setwd("/path/to/your/geode/directory")
   ```

2. **Load the data**:
   ```r
   # Load the preprocessed data
   load("preprocessed_data.RData")
   ```

3. **Run the analysis script**:
   ```r

   # The analysis code is in data_analysis_final.R
   # You can source it or copy relevant sections as needed
   ```

### Key Analysis Sections

The analysis code is organized into several main sections that generate the complete set of figures for a research paper:

#### 1. Data Preprocessing and Variable Creation
- Creates binary group classifications (`binary_group`, `group_verbal`)
- Implements inclusion/exclusion criteria (`srs_normal`, `above_crit`)
- Generates summary statistics and participant tables
- Handles missing data and data quality checks

#### 2. Demographics Analysis (demographics_v2.pdf)
- Comprehensive demographic distributions by group
- Clinical assessment score distributions and group comparisons
- Creates publication-ready demographic plots with proper statistical annotations

#### 3. Main Descriptive Analysis (descriptives_v2.pdf)
- Verbal ability group comparisons and VABS-3 score distributions
- Performance analysis by group with scatter plots
- Learning curve analysis for first 10 trials with exponential fits
- Rolling average analysis for trials 11-200 with confidence intervals

#### 4. Psychometric Curve Analysis (suppl_figure3.pdf)
- Full psychometric curves showing choice behavior vs. flash difference
- Accuracy vs. absolute flash difference analysis
- Learning parameter evolution over time (win-stay, lose-switch)
- Statistical comparisons between groups

#### 5. Temporal Flash Weight Analysis (flash_weights.pdf)
- Flash weight coefficients across 10 temporal bins
- Early, middle, and late temporal window analysis
- Group differences in temporal processing patterns
- Statistical modeling of temporal dynamics

#### 6. Flash Discrimination Heatmaps (heatmaps.pdf, diff_heatmaps.pdf)
- Accuracy heatmaps for TD and ASD groups
- Performance across different flash count combinations
- Group difference visualization
- Numerical discrimination ability assessment

#### 7. Signal Detection Theory Modeling (sdt_model_*.pdf)
- SDT model fitting and parameter estimation
- Sensitivity (d') and bias (c) parameter analysis
- Model fit quality assessment
- Parameter evolution over time
- Linear and non-linear model comparisons

#### 8. Advanced Statistical Modeling (hbm_params.pdf, indglm_params.pdf)
- Hierarchical Bayesian Model parameter estimation
- Individual GLM parameter analysis
- Uncertainty quantification and individual differences

#### 9. Correlation and Dimensionality Analysis (correlations.pdf, cross_correlations.pdf, pca.pdf)
- Inter-variable correlation matrices
- Temporal correlation patterns and lag analysis
- Principal Component Analysis for dimensionality reduction

#### 10. Machine Learning Classification (classify.pdf)
- Feature-based classification using multiple algorithms
- ROC curve analysis for different feature sets
- Classification accuracy metrics and feature importance
- Cross-validation results

#### 11. Test-Retest Reliability (multi_session.pdf)
- Performance consistency across multiple sessions
- Response time reliability analysis
- Learning stability assessment

### Additional Learning-Curve Scripts

The repository also includes focused learning-curve scripts used for bootstrap and subgroup analyses:

- **`scripts/run_learningcurve_bootstrap.R`**  
  Bootstrapped learning-curve analysis for ASD vs TD on the hard phase (trial > 10), fitting  
  `y = A0 + (P − A0) * (1 − exp(L * trial))`.  
  Outputs group-level parameter summaries (A0, P, L), KS tests, bootstrap difference CIs/p-values, and permutation tests.

- **`scripts/run_learningcurve_3group_ks.R`**  
  Three-group bootstrap learning-curve analysis (ASD_high, ASD_low, TD) using the same model and hard-phase trials.  
  Outputs parameter summaries and pairwise KS tests for A0, P, and L across all group pairs.

- **`learningcurveASDTDfit.R`**  
  Earlier exploratory script for ASD vs TD learning-curve fitting and bootstrap comparisons.

- **`scripts/run_nlme_learningcurve.R`**  
  Nonlinear mixed-effects learning-curve models (p0, alpha) by group, fit separately for trials 1–10 and 11–200.  
  Outputs fixed-effect tests and random-effects variance components.

- **`scripts/run_heterogeneity_deep2.R`**  
  Extended heterogeneity analysis: mixed-effects random-slope extraction, variance tests, mixture modeling,  
  and nonlinear mixed-effects fits for trials 1–10 and 11–200.

## Data Quality and Inclusion Criteria

### Primary Inclusion Criteria
- **`srs_normal == 1`**: SRS-2 scores meeting study criteria
- **`game_version == "ft"`**: First-time task completion
- **`above_crit == 1`**: Performance above 61.6% and 200 trials completed

### Group Classifications
- **TD**: Typically developing participants
- **ASD**: Autism spectrum disorder participants
- **ASD (low verbal)**: ASD participants with VABS-3 < 75

### Data Cleaning Notes
- Some participants have multiple sessions (indicated by "_2", "_3" suffixes)
- Missing data is handled explicitly in the analysis
- Performance criteria are applied to ensure data quality

## Output Files and Figures

The analysis code generates a comprehensive set of publication-ready figures for a research paper. All figures are saved in the `results/figures/geodems/` directory:

### Main Figures
- **`demographics_v2.pdf`**: Comprehensive demographic analysis (11 panels: A-K)
  - Gender, age, race, and ethnicity distributions by group
  - Clinical assessment score distributions (SRS-2, BIS, VABS-3, AASP subscales)
  - Publication-ready layout with proper labels and formatting

- **`descriptives_v2.pdf`**: Main descriptive analysis (6 panels: A-F)
  - VABS-3 score distributions by verbal ability group
  - Performance vs. VABS-3 scatter plots
  - Percentage of participants meeting criteria by group
  - Learning curves for first 10 trials
  - Rolling average performance for trials 11-200

### Supplementary Figures
- **`suppl_figure3.pdf`**: Psychometric curve analysis (6 panels: A-F)
  - Full psychometric curves showing choice behavior vs. flash difference
  - Accuracy vs. absolute flash difference
  - Learning parameter evolution over time
  - Win-stay and lose-switch behavior analysis

- **`flash_weights.pdf`**: Temporal flash weight analysis (2 panels: B-C)
  - Flash weight coefficients across 10 temporal bins
  - Early, middle, and late temporal window analysis
  - Group differences in temporal processing

- **`heatmaps.pdf`**: Flash discrimination heatmaps (2 panels: D-E)
  - Accuracy heatmaps for TD and ASD groups
  - Performance across different flash count combinations
  - Visual representation of numerical discrimination abilities

- **`diff_heatmaps.pdf`**: Group difference heatmap (1 panel: F)
  - TD vs. ASD performance differences
  - Areas of relative strength and weakness

### Model-Based Analysis
- **`sdt_model_fit.pdf`**: Signal Detection Theory model fits (2 panels: E-F)
  - Model fit quality assessment
  - Parameter recovery analysis

- **`sdt_model_params.pdf`**: SDT parameter analysis (5 panels: G-K)
  - Sensitivity (d') and bias (c) parameters
  - Parameter evolution over time
  - Group comparisons of model parameters

- **`sdt_model_fit_linear.pdf`**: Linear SDT model analysis
  - Alternative model fitting approaches
  - Model comparison results

- **`sdt_model_params_linear.pdf`**: Linear SDT parameters (5 panels)
  - Linear model parameter estimates
  - Temporal dynamics of parameters

### Advanced Statistical Models
- **`hbm_params.pdf`**: Hierarchical Bayesian Model parameters
  - Individual and group-level parameter estimates
  - Uncertainty quantification

- **`indglm_params.pdf`**: Individual GLM parameter analysis
  - Subject-specific model parameters
  - Individual differences analysis

### Correlation and Dimensionality Analysis
- **`correlations.pdf`**: Correlation matrix analysis
  - Inter-variable correlations
  - Clinical-behavioral relationships

- **`cross_correlations.pdf`**: Cross-correlation analysis
  - Temporal correlation patterns
  - Lag-based relationships

- **`pca.pdf`**: Principal Component Analysis
  - Dimensionality reduction results
  - Component loadings and scores

### Classification Analysis
- **`classify.pdf`**: Machine learning classification results
  - ROC curves for different feature sets
  - Classification accuracy metrics
  - Feature importance analysis

### Test-Retest Reliability
- **`multi_session.pdf`**: Multiple session analysis
  - Test-retest reliability plots
  - Performance consistency across sessions
  - Response time reliability

### Figure Organization
All figures are designed for publication with:
- Consistent color schemes (blue for TD, orange for ASD, red for ASD low-verbal)
- Proper statistical annotations
- Publication-quality formatting
- Multi-panel layouts with clear labels (A, B, C, etc.)
- Standardized themes and styling

## Usage Examples

### Basic Data Exploration
```r
# Load data
load("preprocessed_data.RData")

# View data structure
str(alldata)
dim(alldata)

# Summary statistics
summary(alldata)

# Check group distributions
table(alldata$binary_group)
```

### Subset Data for Analysis
```r
# Get main analysis dataset
main_data <- alldata[game_version == "ft" & srs_normal == 1 & above_crit == 1]

# Get summary by participant
participant_summary <- unique(alldata[, .(username, subID, binary_group, perf, age, gender)])
```

### Create Custom Analyses
```r
# Example: Performance by group
library(ggplot2)
ggplot(main_data, aes(x = binary_group, y = perf)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Group", y = "Performance")

# Example: Flash discrimination analysis
flash_analysis <- alldata[game_version == "ft" & srs_normal == 1 & above_crit == 1, 
                         .(mean_accuracy = mean(score), 
                           mean_rt = mean(rt, na.rm = TRUE),
                           n_trials = .N), 
                         by = .(binary_group, abs(dflash))]

# Example: Temporal flash weight analysis
temporal_weights <- alldata[game_version == "ft" & srs_normal == 1 & above_crit == 1,
                           .(l_early = mean(l_early, na.rm = TRUE),
                             l_middle = mean(l_middle, na.rm = TRUE),
                             l_late = mean(l_late, na.rm = TRUE),
                             r_early = mean(r_early, na.rm = TRUE),
                             r_middle = mean(r_middle, na.rm = TRUE),
                             r_late = mean(r_late, na.rm = TRUE)),
                           by = .(binary_group, group_verbal)]

# Example: Learning strategy analysis (win-stay, lose-switch)
strategy_analysis <- alldata[game_version == "ft" & srs_normal == 1 & above_crit == 1,
                            .(win_stay_rate = mean(winstay, na.rm = TRUE),
                              lose_switch_rate = mean(loseswitch, na.rm = TRUE)),
                            by = .(binary_group, group_verbal)]
```

## Citation and Contact

If you use this dataset in your research, please cite the original publication:
Altered perceptual integration and learning in autism revealed by games inspired by rodent operant tasks
Sucheta Chakravarty, Quan Do, Yutong Li, Helen Tager-Flusberg, Joseph T. McGuire, Benjamin B. Scott
bioRxiv 2025.05.06.652509; doi: https://doi.org/10.1101/2025.05.06.652509

## Troubleshooting

### Common Issues
1. **Package installation errors**: Ensure all required packages are installed
2. **Path issues**: Update the working directory path in the analysis script
3. **Memory issues**: The dataset is ~7.5MB; ensure sufficient RAM for analysis
4. **Missing output directories**: Create `results/figures/geodems/` directory before running

### Data Validation
- Check that `alldata` loads correctly
- Verify group distributions match expected values
- Confirm inclusion criteria are working as intended
