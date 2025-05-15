# # Using Chord Diagram
# install.packages("circlize")
# install.packages("scales")
# install.packages("reshape2")


#load packages
library(circlize)
library(scales)
library(colorspace)
library(reshape2)


#loads in adjacency matrix of MTD values
m <- read.csv("Z:/PRJ-shine_hpc/Aphantasia_fMRI/Figures/Chord Plot/mtd_delta_imagine_all_temp2rsn.csv", header = FALSE, sep = ",")

node_names <- read.csv("Z:/PRJ-shine_hpc/Aphantasia_fMRI/Figures/Chord Plot/labels.csv", header = FALSE, stringsAsFactors = FALSE) #node names
rownames(m) <- node_names[3:16]
colnames(m) <- node_names[1:2]

#convert adjacency matrix into a data frame
m <- t(m)
m <- melt(m)

#create two separate data frames for positive and negative MTD values
mpos <- m
mpos$value[mpos$value < 0] <- 0
mneg <- m
mneg$value[mneg$value > 0] <- 0

mboth <- rbind(mpos,mneg)


# colors <- c(
#   L-TL = "#4d8584",R-TL = "#4d8584",
#  L-M = "#a62f03", R-M = "#a62f03",
#   Rsp = "#00a79c67",
#   Para_Hippo = "#5ac0b9d3",
#   Dorsal_pFC = "#5b5da9",
#   Medial_pFC = "#9799d8",
#   Temp = "#5b5ea9a4",
#   pCun_PCC = "#0681ba"
# )


# chordDiagram(mboth, col = ifelse(mboth$value > 0, "red", "green", order = c("R-TL","L-TL","L-M", "L-FP","L-FT","L-DMN","L-DAN","L-VAN","L-V", "R-M", "R-FP","R-FT","R-DMN","R-DAN","R-VAN","R-V"), big.gap = 10, small.gap = 5))
#svg("Z:/PRJ-shine_hpc/Aphantasia_fMRI/Figures/Chord Plot/chord_diagram.svg", width = 20, height = 20)


chordDiagram(
  mboth, 
  col = ifelse(mboth$value > 0, "#39FF14", "#4B0082"),
  grid.col = "grey",
  annotationTrack = c("grid"),
  preAllocateTracks = 1,
 order = c("R-TL","L-TL","L-M", "L-FP","L-FT","L-DMN","L-DAN","L-VAN","L-V", "R-M", "R-FP","R-FT","R-DMN","R-DAN","R-VAN","R-V"),
  big.gap = 10, small.gap = 5
  ) # plots the diagram as is

# for(si in get.all.sector.index()) {
#   circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2)
# }

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  xplot = get.cell.meta.data("xplot")
  circos.text(mean(xlim), ylim[1], sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(-0, 0))
}, bg.border = NA)