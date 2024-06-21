### I ran it on an interactive node
### srun --partition=phillips --nodes=1 --ntasks-per-node=1 --time=02:00:00 --mem=100000 --account=phillipslab --x11 --pty bash -i
###conda activate genespace
###then in R

library(GENESPACE)
library(ggplot2)


setwd("/gpfs/projects/phillipslab/ateterina/Cbren/phylogenomics/worm23/genespace_newref")


genomes2run <- c("Cbren", "Cbrig","Celeg", "Cinop", "Cnigo","Crema","Ctrop")
gpar <- init_genespace(genomeIDs =c("Cbren", "Cbrig","Celeg", "Cinop", "Cnigo","Crema","Ctrop"), wd="/gpfs/projects/phillipslab/ateterina/Cbren/phylogenomics/worm23/genespace_newref",
path2mcscanx = "/projects/phillipslab/ateterina/scripts/MCScanX")

out <- run_genespace(gpar)



customPal <- colorRampPalette( c("#20567a", "#FF7F00", "#919090", "#1f9186", "#6f6aa1", "#73400e"))
ggthemes <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white"),axis.text.y = element_text(face="italic",size = 10))

ripDat <- plot_riparian(
  gsParam = out,
  palette = customPal,
  refGenome = "Cbren",
  genomeIDs = c("Cinop","Celeg","Ctrop","Cbren","Cbrig","Cnigo","Crema"),
  chrFill = "lightgrey",
  chrExpand = 0.1,
  addThemes = ggthemes,
  forceRecalcBlocks = FALSE,
  labelTheseGenomes = c("Cbren"),chrLabFun = function(x) gsub("^0", "", gsub("chr|scaf|chromosome|scaffold|^lg|_", "", toupper(x)))
  )

p1 <- ripDat$plotData$ggplotObj
p2 <- p1 + scale_y_discrete(limits = c("C. inopinata","C. elegans","C. tropicalis","C. brenneri","C. briggsae","C. nigoni","C. remanei")) + labs(title ="", y ="") + theme(plot.title = element_text(hjust = 0, face = "bold"))

tiff('Fig3.tiff', width = 7.2,height = 5.8, res=300, units = "in") ### Figure 4
p2
dev.off()
