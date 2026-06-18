# ============================================================
# Stratigraphic plotting workshop example
# Site: Prášilské jezero
#
# This script:
# 1. Downloads pollen and diatom dataset metadata from Neotoma
# 2. Extracts pollen/spore sample data
# 3. Calculates relative pollen abundances (%)
# 4. Reshapes the data for stratigraphic plotting
# 5. Produces several rioja/riojaPlot stratigraphic diagrams
# 6. Adds CONISS clustering and DCA ordination curves
# ============================================================

# ============================================================
# 1. Install packages if needed
# ============================================================

# riojaPlot is not on CRAN at the time of writing.
# Install ONE of the following if riojaPlot is not already installed:
# remotes::install_github("nsj3/riojaPlot", build_vignettes = TRUE, dependencies = TRUE)
# install.packages("riojaPlot", repos = "https://nsj3.r-universe.dev")

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
# - keep records with units == "NISP"
#

pra_samp_pollen <- samples(pra_records) %>%
  dplyr::filter(elementtype %in% c("pollen", "spore")) %>%
  dplyr::filter(ecologicalgroup %in% c("TRSH", "UPHE", "VACR")) %>%
  dplyr::filter(units == "NISP")

# ============================================================
# 6. Calculate relative abundance percentages
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
# 7. Create a taxon lookup table
# ============================================================
# This preserves the taxon order by ecological group, then taxon name.
taxa_meta <- pra_pollen_perc %>%
  dplyr::distinct(variablename, ecologicalgroup) %>%
  dplyr::arrange(ecologicalgroup, variablename)


# ============================================================
# 8. Convert long data to wide format
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
# 9. Export/import the dataset to CSV
# ============================================================

data_out<-as.data.frame(t(pra_pollen_wide_perc))

write.table(data_out,
            "pra_pollen.csv",
            sep = ",",
            row.names = TRUE,
            col.names = FALSE)

data_in <- read.table(
  "pra_pollen.csv",
  sep = ",",
  header = FALSE,
  stringsAsFactors = FALSE
)
pra_data <- as.data.frame(t(data_in))
colnames(pra_data) <- pra_data[1, ]
pra_data <- pra_data[-1, ]
rownames(pra_data) <- seq_len(nrow(pra_data))
pra_data <- lapply(pra_data, as.numeric)

# ============================================================
# 10. Prepare plotting objects
# ============================================================

# Taxa matrix
taxa <- pra_pollen_wide_perc %>%
  dplyr::select(-depth, -age) %>%
  as.matrix()

# Sample age and depth vectors
age <- pra_pollen_wide_perc$age
depth <- pra_pollen_wide_perc$depth

# ============================================================
# 11. Remove very rare taxa for readability
# ============================================================
# Because there are many taxa, the diagram can become overcrowded.
# If you want a stricter or looser filter, change the threshold.

mx <- apply(taxa, 2, max, na.rm = TRUE)
spec <- taxa[, mx > 2, drop = FALSE]

# ============================================================
# 12. Basic stratigraphic plots
# ============================================================
# Basic stratigraphic plot
strat.plot(
  spec
    )

# Percentage-scaled plot
strat.plot(
  spec,
  y.rev = TRUE,
  scale.percent = TRUE,
  ylim = c(0, 200) # Y-axis limits
)

# Plot by sample depth and inpercentage scale
strat.plot(
  spec,
  yvar = depth,
  y.rev = TRUE,
  ylim = c(1480, 1700),
  ylabel = "Depth (cm)",
  scale.percent = TRUE
)

# Filled polygons and sample bars
strat.plot(
  spec,
  yvar = depth,
  y.rev = TRUE,
  ylim = c(1480, 1700),
  ylabel = "Depth (cm)",
  scale.percent = TRUE,
  plot.poly = TRUE,
  col.poly = "darkgreen",
)

# Add exaggeration curves
strat.plot(
  spec,
  yvar = depth,
  y.rev = TRUE,
  ylabel = "Depth (cm)",
  scale.percent = TRUE,
  plot.poly = TRUE,
  col.poly = "darkgreen",
  exag = TRUE,
  exag.mult = 5,
  col.exag = "yellow"
)


# ============================================================
# 13. Colour curves by taxon type
# ============================================================

# This is a simple manual colouring example.
# Adjust the break points depending on how your taxa are ordered.
xx <- 1:ncol(spec)
cc <- ifelse(xx < 14, "darkgreen",
             ifelse(xx > 20, "blue", "red"))

# Plot by calibrated age
strat.plot(
  spec,
  yvar = age,
  y.rev = TRUE,
  y.tks = c(0, 1000, 2000, 3000, 4000, 5000, 6000,
            7000, 8000, 9000, 10000, 11000),
  ylabel = "Age (cal BP)",
  ylabPos = 3,
  scale.percent = TRUE,
  cex.yaxis = 0.8,
  tcl = -0.25,
  plot.poly = TRUE,
  col.poly = cc,
  exag = TRUE,
  exag.mult = 5,
  col.exag = "yellow"
)

# ============================================================
# 14. Plot taxon names in italic
# ============================================================
# 1. Define the exact names you want to KEEP AS PLAIN TEXT
exclude_names <- c("Cyperaceae", "Poacese", "Poaceae (Cerealia)", "Polypodiaceae")

# 2. Convert names conditionally
italic_labels <- as.expression(lapply(colnames(spec), function(x) {
  if (x %in% exclude_names) {
    return(x) # Returns plain text (no italics)
  } else {
    return(bquote(italic(.(x)))) # Returns italic expression
  }
}))

# 3. Run your exact base strat.plot code, substituting the names
strat.plot(
  spec,
  x.names = italic_labels, # This is what turns them italic!
  yvar = depth,
  y.rev = TRUE,
  ylabel = "Depth (cm)",
  scale.percent = TRUE
)

# ============================================================
# 15. Constrained cluster analysis (CONISS)
# ============================================================

# Transform the dataset and calculate dissimilarity and run stratigraphically constrained clustering.

diss <- dist(sqrt(taxa))
clust <- chclust(diss, method = "coniss")

# Broken-stick model to help decide the number of significant zones
bstick(clust)

# Plot with dendrogram
x<-strat.plot(
  spec,
  x.names = italic_labels, # This is what turns them italic!
  min.width=10,
  yvar = age,
  y.rev = TRUE,
  y.tks = c(0, 1000, 2000, 3000, 4000, 5000, 6000,
            7000, 8000, 9000, 10000, 11000),
  scale.percent = TRUE,
  cex.yaxis = 0.8,
  ylabel = "Age (cal BP)",
  ylabPos = 3,
  plot.poly = TRUE,
  col.poly = cc,
  exag = TRUE,
  exag.mult = 5,
  col.exag = "yellow",
  clust = clust
)

# Add zone boundaries
# Here, "6" means draw the boundary corresponding to 6 cuts / a chosen level.
# Adjust depending on the number of zones suggested by broken stick model.
addClustZone(x, clust, 6, col = "red")

# ============================================================
# 16. Combine pollen diagram with DCA sample scores
# ============================================================

# Use fig/xLeft/xRight controls to place multiple diagrams in one figure.

strat.plot(
  spec,
  xRight = 0.7,
  yvar = depth,
  y.rev = TRUE,
  scale.percent = TRUE,
  ylabel = "Depth (cm)"
)

# Detrended Correspondence Analysis (DCA) on pollen percentages
dca <- decorana(spec, iweigh = 1)

# Extract site/sample scores for axes 1 and 2
sc <- scores(dca, display = "sites", choices = 1:2)

# Add DCA curves to the right side of the pollen diagram
strat.plot(
  sc,
  xLeft = 0.7,
  xRight = 1,
  yvar = depth,
  y.rev = TRUE,
  y.axis = FALSE,
  clust = clust,
  clust.width = 0.08,
  add = TRUE
)

# ============================================================
# 17. Add a custom smoothing function
# ============================================================

# This function adds a LOWESS smooth to each taxon curve.
# In strat.plot(), custom functions receive:
# - x: taxon values
# - y: stratigraphic variable (e.g. depth)
# - i: taxon index
# - nm: taxon name

sm.fun <- function(x, y, i, nm) {
  tmp <- data.frame(x = y, y = x)
  tmp <- na.omit(tmp)
  
  lo <- lowess(tmp, f = 0.1)
  lines(lo$y, lo$x, col = "red", lwd = 1)
}

# Plot using the custom smoothing function
strat.plot(
  spec,
  yvar = depth,
  y.rev = TRUE,
  ylabel = "Depth (cm)",
  scale.percent = TRUE,
  fun1 = sm.fun
)

# ============================================================
# 18. Save the diagram as PDF 
# ============================================================

pdf("strat_plot.pdf",
    width = 8,
    height = 10)

strat.plot(
  spec,
  yvar = depth,
  y.rev = TRUE,
  ylabel = "Depth (cm)",
  scale.percent = TRUE,
  fun1 = sm.fun
)

dev.off()
