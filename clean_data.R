# Read in R packages
# Note: if you don't have some of these installed, you can do so by running
# install.packages("package_name")

library(readr)
library(dplyr)
library(openxlsx)

# ============================================
# Data Cleaning

# We need to roll up rows that share key identifiers: species name, edibility,
# etc. When we roll those up, we can simply sum the values for columns J to N.


# ============================================
# Write Files to Disk


