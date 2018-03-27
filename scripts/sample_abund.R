# *************************************************************************
# Script to sample abundance considering a minimum area factor.
# *************************************************************************

library(data.table)  # for data manipulation & aggregation
library(magrittr)    # for Unix-like pipeline syntax
library(matrixStats) # for matrix fast operations

# Load `do_sample()` function
source("scripts/do_sample.R")


# Read & prepare data -----------------------------------------------------

dt <- fread("data/Hawaii_Plot_SppAbundance_Filter.csv")

# Select only the first row from each plot.
dt_unq <- unique(dt, by = "PlotIDn")


# Data exploration --------------------------------------------------------

# Histogram of plot area
dt_unq[, Plot_Area] %>% hist(breaks = 200)
# Histogram of plot area only for plots smaller than ... 
dt_unq[Plot_Area < 1200, Plot_Area] %>% hist(breaks = 200)

# Counts of plots per islands, for plots of given area constraints
dt_unq[Plot_Area > 1000, .N, by = Island]
dt_unq[Plot_Area > 500, .N, by = Island]
dt_unq[Plot_Area %between% c(400,600), .N, by = Island]
dt_unq[Plot_Area %between% c(600,800), .N, by = Island]
dt_unq[Plot_Area %between% c(900, 1200), .N, by = Island]

rm(dt_unq) # no further need of dt_unq


# Sample ------------------------------------------------------------------

# STEPS

# Subset data as needed
# Get the total abundance by plot
# Explode rows using raw abundance
# Group exploded rows by plot ID and:
# - compute specific sample size considering area factor
# - sample n times the species names with the specific sample size
# - for each sample, get counts by species, so get a counts matrix
# - from the counts matrix, get average of counts and SD


# Use a subset of data for sampling
dt2 <- dt[Plot_Area > 600, .(PlotIDn, Abundance_raw, SppCode, Plot_Area)]
min_area <- min(dt2$Plot_Area) # min sample area
# dt2[, min(Plot_Area)] # or use canonical data.table syntax

# Get the total abundance by plot
dt2[, sum_abund := sum(Abundance_raw), by = PlotIDn]

# Explode rows using raw abundance
dt2_expl <- dt2[rep(1:.N, times = Abundance_raw)]
# the line below is valid both for data.frame and data.table
# dt2_expl <- dt2[rep(seq_len(nrow(dt2)), times = dt2$Abundance_raw),]

# Open the exploded data table and, while grouping by plot ID, do the sampling
sample_results <- dt2_expl[, do_sample(sp_code = SppCode, 
                                       min_area = min_area,
                                       plot_area = Plot_Area, 
                                       sum_abund = sum_abund, 
                                       repl = 1000),
                           by = PlotIDn]
# ~ 40 sec for 1000 replicates

results <- merge(x = dt[Plot_Area > 600],
                 y = sample_results,
                 by.x = c("PlotIDn", "SppCode"),
                 by.y = c("PlotIDn", "sp_code"),
                 sort = FALSE)
View(results)

write.csv(results, 
          file = "output/Hawaii_Plot_SppAbundance_min_area_sampled.csv", 
          row.names = FALSE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Optional - `for` loop ---------------------------------------------------

# The 'sample_results' object can also be constructed with the traditional 'for' loop aproach.
# Is a bit more verbose but maybe could be easier to fallow.
# Also it may be easier to debug when the code fails at a certain step (i)

# Split the exploded table by plot id
dt_lst <- split(dt2_expl, by = "PlotIDn")

# Allocate memory for an empty list
results_lst <- vector(mode = "list", length = length(dt_lst))

# For each plot execute the sampling procedure
for (i in 1:length(dt_lst)){
  spl <- do_sample(sp_code = dt_lst[[i]]$SppCode, 
                   min_area = min_area,
                   plot_area = dt_lst[[i]]$Plot_Area, 
                   sum_abund = dt_lst[[i]]$sum_abund, 
                   repl = 1000)
  results_lst[[i]] <- as.data.table(spl)
}

sample_results <- rbindlist(results_lst)
# then do the merging
