# EigenRF: An Improved Metabolomics Normalization Method

# Overview
EigenRF is an R package designed to improve the normalization of metabolomics data. It enhances the previous EigenMS method by incorporating a random forest regression model to capture nonlinear biological variations of interest. This method not only eliminates systematic errors but also preserves the biological variations of interest, leading to improved accuracy and reproducibility in metabolomics research.

# Installation
You can install the EigenRF package directly from GitHub using the following command:

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("statisticsinMS/EigenRF")

# Usage
To use the EigenRF normalization method, simply load the library and apply the Eigen_RF function to your data:

library(EigenRF)

# Your peak table (e.g., a data frame with peak intensities)
peak <- your_peak_data

# Your group vector (e.g., a factor indicating the group of each sample)
groups <- your_group_data

# Your metabolites table or vector
metabolites <- your_metabolites_data

# Apply the Eigen_RF normalization
normalized_data <- Eigen_RF(peak, groups, metabolites)

# Function Parameters
peak: A data frame containing the peak intensity values.
groups: A vector indicating the group affiliation of each sample.
metabolites: A data frame or vector containing metabolite data.

# Details
The Eigen_RF function returns a list with the normalized peak table, facilitating the analysis of differential metabolites with improved reproducibility.

# Acknowledgements
This research was funded by the National Key R&D Program of China and the National Natural Science Foundation of China.

# Authors:
Chencheng Tang, Dongfang Huang, Xudong Xing, Hua Yang.
