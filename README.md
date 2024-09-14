
# HLA Gene Analysis

This project involves the extraction, analysis, and classification of HLA-A and HLA-B gene sequences using bioinformatics and machine learning techniques. The analysis includes data acquisition from NCBI, sequence filtering, nucleotide content analysis, and classification using Random Forest and Naive Bayes algorithms.

## Repository Structure

```
HLA_Gene_Analysis/
│
├── data/                        # Directory for input data files
│   ├── HLA_A_seq.fasta
│   └── HLA_B_seq.fasta
│
├── scripts/                     # Directory for R scripts
│   └── HLA_analysis.R
│
├── output/                      # Directory for output files like plots and reports
│
├── README.md                    # README file
├── .gitignore                   # Gitignore file
└── LICENSE                      # License file
```

## Requirements

- R (version >= 4.0)
- R Packages:
  - BiocManager
  - Biostrings
  - bioseq
  - qpcR
  - rgl
  - glmnet
  - e1071
  - tidyverse
  - ggplot2
  - reshape2
  - data.table
  - randomForest

## Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/HLA_Gene_Analysis.git
   cd HLA_Gene_Analysis
   ```

2. **Install R and R packages**:

   Open R or RStudio and run the following script to install the required packages:
   ```r
   # Install Bioconductor and other required packages
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install("Biostrings", version = "3.15", force = TRUE)
   install.packages(c('bioseq', 'qpcR', 'rgl', 'glmnet', 'e1071', 'tidyverse', 'ggplot2', 'reshape2', 'data.table', 'randomForest'))
   ```

3. **Run the Analysis**:

   Navigate to the `scripts` directory and run the R script:
   ```bash
   Rscript scripts/HLA_analysis.R
   ```

   This script will perform the analysis and save the output files in the `output` directory.

## Analysis Overview

- **Data Acquisition**: Fetch HLA-A and HLA-B sequences from the NCBI database.
- **Data Filtering**: Remove incomplete and unusually short sequences.
- **Sequence Analysis**: Calculate nucleotide frequencies and visualize GC and AT content.
- **Classification**: Use Random Forest and Naive Bayes classifiers to categorize sequences.
- **Visualization**: Generate plots to illustrate the nucleotide content and classifier performance.

## Output

- `output/Frequency_Boxplot.png`: Boxplot of AT & GC frequency.
- `output/Mean_length_barplot.png`: Barplot of the average length of HLA-A and HLA-B sequences.
- `output/Random_Forest_Plot.png`: Confusion matrix for Random Forest classification.
- `output/Naive_Bayes_Plot.png`: Confusion matrix for Naive Bayes classification.
