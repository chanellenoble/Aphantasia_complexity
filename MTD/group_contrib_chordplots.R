# Citation: Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize Implements and enhances circular visualization in R. Bioinformatics. 2014 Oct;30(19):2811-2. doi: 10.1093/bioinformatics/btu393. Epub 2014 Jun 14. PMID: 24930139

library(circlize)
library(scales)
library(reshape2)
library(dplyr)

# Function to load edges CSV and plot chord diagram
plot_chord <- function(filename, title_text) {
  edges <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
  colnames(edges) <- c("from", "to", "value")
  
  # Clean 'to' labels: remove trailing underscore+number to map to RSN only
  edges$to <- gsub("_\\d+$", "", edges$to)
  
  # Map RSN codes to full names
  rsn_rename <- c(
    RSN1 = "L-V",   RSN2 = "L-M",   RSN3 = "L-DAN", RSN4 = "L-VAN",
    RSN5 = "L-FT",  RSN6 = "L-FP",  RSN7 = "L-DMN", RSN8 = "R-V",
    RSN9 = "R-M",   RSN10 = "R-DAN", RSN11 = "R-VAN", RSN12 = "R-FT",
    RSN13 = "R-FP", RSN14 = "R-DMN"
  )
  
  # Fix TL labels and apply RSN name mapping
  edges$from <- ifelse(edges$from == "LT", "L-TL",
                       ifelse(edges$from == "RT", "R-TL", edges$from))
  edges$to <- ifelse(edges$to %in% names(rsn_rename),
                     rsn_rename[edges$to],
                     edges$to)
  
  # Define your preferred sector order
  sector_order <- c("R-TL","L-TL","L-M", "L-FP","L-FT","L-DMN","L-DAN","L-VAN","L-V","R-V","R-VAN", "R-DAN","R-DMN","R-FT" , "R-FP", "R-M")
  
  # Clean whitespace to ensure exact string matching
  edges$from <- trimws(as.character(edges$from))
  edges$to <- trimws(as.character(edges$to))
  sector_order <- trimws(as.character(sector_order))
  
  # Drop edges not in sector_order to avoid circos error
  edges <- edges[edges$from %in% sector_order & edges$to %in% sector_order, ]
  
  # Check for any mismatches
  all_labels <- unique(c(edges$from, edges$to))
  missing_sectors <- setdiff(all_labels, sector_order)
  if (length(missing_sectors) > 0) {
    warning("Missing sectors in order: ", paste(missing_sectors, collapse = ", "))
  }
  
  # Assign edge colors
  edges$col <- ifelse(edges$value > 0, "#39FF14", "#4B0082")
  
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
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        x = mean(xlim),
        y = ylim[1],
        labels = sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.2
      )
    },
    bg.border = NA
  )
  
  title(main = title_text, cex.main = 1.5)
}

# Run for each group
plot_chord("C:/Users/Chanelle/OneDrive - The University of Sydney (Students)/Sydney Uni Work/Neuroscience/Shine-Labs/Aphantasia Manuscript/Figures-2025-07-24/Circlize/edges_for_circlize_aphan.csv", "Aphantasics - MTD Beta Contributions")
plot_chord("C:/Users/Chanelle/OneDrive - The University of Sydney (Students)/Sydney Uni Work/Neuroscience/Shine-Labs/Aphantasia Manuscript/Figures-2025-07-24/Circlize/edges_for_circlize_ctrl.csv",  "Controls - MTD Beta Contributions")
plot_chord("C:/Users/Chanelle/OneDrive - The University of Sydney (Students)/Sydney Uni Work/Neuroscience/Shine-Labs/Aphantasia Manuscript/Figures-2025-07-24/Circlize/edges_for_circlize_hyper.csv", "Hyperphantasics - MTD Beta Contributions")

