# ============================================================
# Stratigraphic plotting workshop example
# Site: Prášilské jezero
#
# This script:
# 1. Downloads pollen and diatom dataset metadata from Neotoma
# 2. Extracts pollen/spore sample data
# 3. Calculates relative pollen abundances (%)
# 4. Reshapes the data for stratigraphic plotting
# 5. Produces  rioja/riojaPlot stratigraphic diagrams
# 6. Adds CONISS clustering
# ============================================================


# ============================================================
# 1. Install packages if needed
# ============================================================

# riojaPlot is not on CRAN at the time of writing.
# Install ONE of the following if riojaPlot is not already installed:

# remotes::install_github("nsj3/riojaPlot", build_vignettes = TRUE, dependencies = TRUE)
# install.packages("riojaPlot", repos = "https://nsj3.r-universe.dev")

# Optional: install other packages if needed
# install.packages(c("tidyverse", "vegan", "rioja", "sjmisc"))


# ============================================================
# 2. Load libraries
# ============================================================

library(neotoma2)
library(tidyverse)
library(rioja)
library(riojaPlot)
library(vegan)

# ============================================================
# 3. Get dataset metadata for the site
# ============================================================

pra_site <- get_datasets(get_sites(sitename = "Prášilské jezero"), all_data = TRUE) %>%
  neotoma2::filter(datasettype %in% c("pollen","diatom"))


# ============================================================
# 4. Download records and inspect contents
# ============================================================
# Download the selected datasets
pra_records <- pra_site %>% get_downloads()

# View the datasets object
datasets(pra_records)

# Extract sample-level information
pra_samp <- pra_records %>% samples()

# Extract taxa metadata
pra_taxa <- neotoma2::taxa(pra_records)


# ============================================================
# 5. Filter pollen/spore data for pollen sum
# ============================================================
# We want upland pollen and spores for the pollen sum.
# Ecological groups used here:
# - TRSH = trees and shrubs
# - UPHE = upland herbs
# - VACR = vascular cryptogams / ferns
#
# We also:
# - keep only pollen and spores
# - exclude stomata by filtering elementtype
# - keep records containg Number of Identified Speciments == "NISP"

pra_samp_pollen <- samples(pra_records) %>%
  dplyr::filter(ecologicalgroup %in% c("TRSH", "UPHE", "VACR")) %>%
  dplyr::filter(elementtype %in% c("pollen", "spore")) %>%
  dplyr::filter(units == "NISP")

# ============================================================
# 6. Store taxon-to-ecological-group lookup
# ============================================================
# This table is useful later if you want to group taxa
# or annotate the diagram by ecological category.
pra_ecol_groups <- pra_samp_pollen %>%
  dplyr::select(variablename, ecologicalgroup) %>%
  distinct()


# ============================================================
# 7. Calculate relative abundance percentages
# ============================================================
# For each sample (defined here by depth and age), calculate the total
# pollen sum and then express each taxon as a percentage of that sum.

pra_pollen_perc <- pra_samp_pollen %>%
  dplyr::group_by(depth, age) %>%
  dplyr::mutate(pollen_sum = sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(relative_abundance_percent = (value / pollen_sum) * 100) %>%
  dplyr::arrange(dplyr::desc(age))

# ============================================================
# 8. Create a taxon lookup table
# ============================================================
# This arranges the taxon order by ecological group, then taxon name.
taxa_meta <- pra_pollen_perc %>%
  dplyr::distinct(variablename, ecologicalgroup) %>%
  dplyr::arrange(ecologicalgroup, variablename)

# ============================================================
# 9. Convert long data to wide format
# ============================================================
# rioja::strat.plot() expects a matrix/data frame with:
# - rows = samples
# - columns = taxa
#
# Here we reshape the percentage data so each taxon becomes a column.
# Missing taxa in a sample are filled with 0.

pra_pollen_wide_perc <- pra_pollen_perc %>%
  dplyr::select(depth, age, variablename, relative_abundance_percent) %>%
  tidyr::pivot_wider(
    id_cols = c(depth, age),
    names_from = variablename,
    values_from = relative_abundance_percent,
    values_fill = 0
  ) %>%
  dplyr::arrange(depth) %>%
  dplyr::select(depth, age, dplyr::all_of(taxa_meta$variablename))


# ============================================================
# 10. Prepare plotting objects - these should be dibbles or data frames!
# ============================================================

# Taxa matrix for plotting
taxa <- pra_pollen_wide_perc %>%
  dplyr::select(-depth, -age) %>%
  as.data.frame()

chron <- pra_pollen_wide_perc %>%
  dplyr::select(depth, age) %>%
  dplyr::rename(
    'Depth (cm)' = depth,
    `Age (cal BP)` = age
  ) %>%
  as.data.frame()


# ============================================================
# Start with riojaplot!
# ============================================================
# ---------------------------------------------------------
# 1. Minimal plot
# riojaPlot takes:
#   x = variables to plot
#   y = chronology/depth/time data
# ---------------------------------------------------------
riojaPlot(taxa, chron)


# ---------------------------------------------------------
# 2. Percentage-scaled plot
# Useful for pollen percentage data
# ---------------------------------------------------------
riojaPlot(taxa, chron, scale.percent = TRUE)

# ---------------------------------------------------------
# 3. Plot only common taxa
# Keep taxa with maximum abundance > 5
# ---------------------------------------------------------
max5 <- sapply(taxa, max) > 2
sel_taxa <- names(max5[max5])

riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  scale.percent = TRUE
)

# ---------------------------------------------------------
# 4. Compare styles
# riojaPlot can show silhouettes, lines, and/or bars
# ---------------------------------------------------------

# Default style
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  scale.percent = TRUE
)

# Silhouettes only
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  scale.percent = TRUE,
  plot.line = FALSE,
  plot.bar = FALSE
)

# Lines only
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  scale.percent = TRUE,
  plot.poly = FALSE,
  plot.bar = FALSE
)

# Bars only
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  scale.percent = TRUE,
  plot.line = FALSE,
  plot.poly = FALSE,
  plot.bar = TRUE,
  lwd.bar = 2
)


# ---------------------------------------------------------
# 5. Add groups
# Use grouping information from aber$names
# ---------------------------------------------------------
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  groups = taxa_meta,
  scale.percent = TRUE,
  plot.groups = TRUE
)

# ---------------------------------------------------------
# 6. Add cumulative group plot
# Summarises groups such as Trees / Shrubs / Herbs
# ---------------------------------------------------------
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  groups = taxa_meta,
  scale.percent = TRUE,
  plot.groups = TRUE,
  plot.cumul = TRUE
)

# ---------------------------------------------------------
# 7. Add exaggeration curves
# Helps show low-abundance taxa more clearly
# ---------------------------------------------------------
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  scale.percent = TRUE,
  plot.exag = TRUE
)


# ---------------------------------------------------------
# 8. Add a secondary y-axis
# Example: depth and age scales together
# ---------------------------------------------------------
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  yvar.name = "Depth (cm)",
  sec.yvar.name = "Age (cal BP)",
  plot.sec.axis = TRUE,
  scale.percent = TRUE,
  plot.exag = TRUE
)

# ---------------------------------------------------------
# 9. Add zonation automatically
# Cluster the data and draw zone boundaries
# ---------------------------------------------------------
riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  groups = taxa_meta,
  scale.percent = TRUE,
  plot.groups = TRUE,
  plot.cumul = TRUE,
  clust.data.trans = "sqrt",
  do.clust = TRUE,
  plot.clust = TRUE,
  plot.zones = "auto"
)

# ---------------------------------------------------------
# 10. Add zones manually
# Save the plot object, then annotate it
# ---------------------------------------------------------
rp <- riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  yvar.name = "Age (cal BP)",
  ymin = 0,
  ymax = 11300,
  yinterval = 500,
  scale.percent = TRUE,
  cex.xaxis = 0.5,
  cex.xlabel = 0.7,
  las.xaxis = 2
)

# Add horizontal zone lines
addRPZone(rp, c(8200, 4200), col = "red")

# Add shaded zones
addRPZone(rp, 5000, 78000)
addRPZone(rp, 1000, 1500, col = "blue", alpha = 0.05)


# ---------------------------------------------------------
# 11. Save figure directly to file
# Drawing straight to svg() usually gives better output
# than resizing from the plot window
# ---------------------------------------------------------
pdf("riojaPlot_workshop_example.pdf", width = 8, height = 5)

riojaPlot(
  taxa,
  chron,
  selVars = sel_taxa,
  yvar.name = "Depth (cm)",
  sec.yvar.name = "Age (cal BP)",
  plot.sec.axis = TRUE,
  scale.percent = TRUE,
  plot.exag = TRUE
)

dev.off()

