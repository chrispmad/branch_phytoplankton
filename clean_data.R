# Read in R packages
# Note: if you don't have some of these installed, you can do so by running
# install.packages("package_name")

library(readr)
library(dplyr)
library(openxlsx)

# ============================================
# Read in Data files

al = readr::read_csv("data/ALU_DAT_2023-2024_phyto-comparison_forCM.csv")
wa = readr::read_csv("data/WAH_DAT_2023-2024_phyto-comparison_2025-04-11_forCM.csv")

# We can drop 'day' from 'wa' - this info is contained in the date column.
wa = wa |>
  dplyr::select(-day)

# Combine data files, adding a column noting which project each row is from.
al$project = 'Alouette'
wa$project = 'Wahleach'

d = dplyr::bind_rows(al, wa) |>
  # Rearrange the order of the columns: put 'project' first, then all the remaining columns.
  dplyr::select(project, dplyr::everything())

# ============================================
# Data Cleaning

# Some very minor typos to clean up:
# 1. class.alias
d = d |>
  dplyr::mutate(class.alias = stringr::str_replace_all(class.alias,";",","))

# We need to roll up rows that share key identifiers: species name, edibility,
# etc. When we roll those up, we can simply sum the values for columns J to N.
d_rolled = d |>
  # First, let's drop any extraneous detail in the species column, e.g. size of organism.
  dplyr::mutate(spp = stringr::str_remove(spp, "\\((small|medium|large)\\)")) |>
  # We're going to group by most of the columns - we don't want to lose that granularity
  # when we roll up the lab results.
  dplyr::group_by(project,date,year,month,basin,type,
                  depth.m,class.alias,class.name,spp,edibility) |>
  # Now we'll apply the rolling-up to all four numeric columns, summing them.
  dplyr::reframe(
    dplyr::across(
      c(DB.cells.ml:LV.mm3.l), \(x) sum(x,na.rm=T)
    )
  )

print(paste0("We've condensed the data from ",nrow(d)," rows to ",nrow(d_rolled), " rows."))

# ============================================
# Write Files to Disk
openxlsx::write.xlsx(x = d_rolled, "output/data_both_projects_rolled_up.xlsx", overwrite = T)

