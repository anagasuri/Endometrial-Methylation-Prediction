# Endometrial Methylation Prediction

## Overview

This project investigates whether DNA methylation profiles from eutopic endometrial tissue can be used to predict endometriosis status and menstrual cycle phase. Using Illumina EPIC methylation array data from 984 endometrial samples, I developed machine learning pipelines to evaluate predictive performance, identify informative CpG sites, and explore the biological pathways associated with disease status and normal endometrial cycling.

This work was completed as part of my M.S. in Health Data Science at UCSF and serves as the basis for my graduate capstone project.

## Research Questions

1. Can DNA methylation profiles accurately distinguish endometriosis cases from controls?
2. How strongly does menstrual cycle phase influence endometrial methylation patterns?
3. Which CpG sites contribute most to prediction performance?
4. Do predictive methylation signatures reveal biologically meaningful pathways?

## Dataset

### Cohort

- 984 eutopic endometrial biopsy samples
- 637 endometriosis cases
  - 344 Stage I/II
  - 286 Stage III/IV
  - 7 unknown stage
- 347 controls
- Samples collected across multiple international cohorts:
  - UCSF
  - Oxford
  - Melbourne
  - Edinburgh

### Data Types

#### DNA Methylation

- Illumina EPIC BeadChip Array
- 759,345 CpG sites after quality control

#### Clinical Metadata

- Age
- BMI
- Parity
- Age at menarche
- Genetic ancestry
- Menstrual cycle phase
- Disease stage
- Lesion type
- Pain phenotypes

## Workflow

### 1. Batch Correction

Technical variation resulting from batch effects and processing site differences was corrected using SmartSVA.

Quality assessment was performed using:

- Principal Component Analysis (PCA)
- Uniform Manifold Approximation and Projection (UMAP)
- Silhouette score analysis

SmartSVA substantially reduced clustering driven by batch and processing site while preserving biologically meaningful variation.

### 2. Full-Feature Modeling

Ridge Logistic Regression (L2 regularization) models were trained using all 759,345 CpG sites.

Two classification tasks were evaluated:

- Endometriosis Case vs. Control
- Menstrual Cycle Phase

Performance was assessed using:

- AUROC
- Precision-Recall Curves
- Held-out test sets
- Label permutation experiments

### 3. Feature Selection

Several feature-selection strategies were explored:

#### Elbow Method

- CpGs ranked by absolute model coefficient magnitude
- Knee-point detection used to identify informative feature subsets

#### T-Test

- Differential methylation testing
- Benjamini-Hochberg FDR correction

#### Univariate Logistic Regression

- Individual CpG association testing
- Multiple-testing correction

### 4. Biological Interpretation

Selected CpGs were evaluated using:

- Gene Ontology (GO)
- KEGG Pathways
- Reactome Pathways
- Heatmap visualization

## Key Results

### Technical Batch Effects Were Successfully Reduced

PCA and UMAP analyses demonstrated substantial reduction in clustering by:

- Batch
- Institute of analysis

following SmartSVA correction.

Silhouette scores confirmed reduced technical structure while preserving biological signal.

### Menstrual Cycle Phase Is the Dominant Source of Methylation Variation

Cycle phase produced substantially stronger methylation signals than disease status throughout both exploratory analyses and predictive modeling.

This finding supports previous observations that normal endometrial biology is a major driver of methylation variation.

### Machine Learning Accurately Predicts Endometriosis Status

Using genome-wide methylation profiles, ridge logistic regression models achieved strong predictive performance for distinguishing endometriosis cases from controls.

Permutation testing demonstrated that predictive performance was not driven by random chance.

### Machine Learning Accurately Predicts Menstrual Cycle Phase

Models achieved even stronger performance when predicting proliferative versus secretory phase samples.

These results indicate that DNA methylation profiles contain highly informative signatures of endometrial cycling.

### Feature Selection Improves Interpretability

Multiple feature-selection strategies were compared to identify CpG subsets associated with prediction performance.

The elbow method provided compact feature sets that maintained predictive accuracy while improving biological interpretability.

### Biological Pathways Reflect Endometrial Function

Pathway enrichment analyses identified numerous biologically relevant pathways associated with cycle-phase-related methylation signatures, while disease-associated signals appeared more diffuse.

## Tools and Technologies

### Programming Languages

- Python
- R

## Repository Contents

- Data preprocessing workflows
- SmartSVA batch correction pipelines
- PCA and UMAP analyses
- Machine learning models
- Feature-selection workflows
- Permutation testing analyses
- Pathway enrichment analyses
- Heatmap generation scripts
- Figure-generation notebooks

## Author

**Amrita Nagasuri**

M.S. Health Data Science  
University of California, San Francisco





