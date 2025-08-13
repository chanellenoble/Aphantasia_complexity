# Load required packages
library(circlize)
library(scales)
library(reshape2)
library(dplyr)
# Load the edge list CSV exported from MATLAB
edges <- read.csv("C:/Users/Chanelle/OneDrive - The University of Sydney (Students)/Sydney Uni Work/Neuroscience/Shine-Labs/Aphantasia Manuscript/Figures-2025-07-24/Circlize/edges_for_circlize.csv")

colnames(edges) <- c("from", "to", "value")

# Clean up RSN labels: collapse "RSN1_3" to just "RSN1"
edges$to <- gsub("_\\d+$", "", edges$to)



# Map RSN codes to full names (customize this list as needed!)
rsn_rename <- c(
  RSN1 = "L-V",
  RSN2 = "L-M",
  RSN3 = "L-DAN",
  RSN4 = "L-VAN",
  RSN5 = "L-FT",
  RSN6 = "L-FP",
  RSN7 = "L-DMN",
  RSN8 = "R-V",
  RSN9 = "R-M",
  RSN10 = "R-DAN",
  RSN11 = "R-VAN",
  RSN12 = "R-FT",
  RSN13 = "R-FP",
  RSN14 = "R-DMN"
)


edges$to <- rsn_rename[edges$to]

# Optional: collapse L and R temporal lobes into one label each
edges$from <- ifelse(edges$from == "LT", "L-TL",
                     ifelse(edges$from == "RT", "R-TL", edges$from))


# Set final order of sectors manually
sector_order <- c("R-TL","L-TL","L-M", "L-FP","L-FT","L-DMN","L-DAN","L-VAN","L-V","R-V","R-VAN", "R-DAN","R-FT" , "R-FP", "R-M")


# Color positive vs negative correlations
edges$col <- ifelse(edges$value > 0, "#39FF14", "#4B0082")

# Define all sectors for consistent ordering
all_sectors <- unique(c(edges$from, edges$to))

edges$from <- factor(edges$from, levels = sector_order)
edges$to <- factor(edges$to, levels = sector_order)

# Plot
circos.clear()
chordDiagram(
  edges[, c("from", "to", "value")],
  col = edges$col,
  grid.col = "grey",
  annotationTrack = "grid",
  preAllocateTracks = 1,
  order = sector_order,
  big.gap = 12,
  small.gap = 1
)

# Add labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2)
}, bg.border = NA)