##############################################
######### C.breneri genome ###################
########### Figures ##########################
##############################################

setwd("/Users/anastasia/Documents/Phillips_lab/drafts/Cbren_genome")

library(ggplot2)

##############################################
################## Fig2 ######################
##############################################

#library(ggplot2)
EXON<-read.csv("Cbren_Exons_100Kb.bed",header=FALSE, sep="\t")
INTRON<-read.csv("Cbren_Introns_100Kb.bed",header=FALSE, sep="\t")
REP<-read.csv("Cbren_Repeats_100Kb.bed",header=FALSE, sep="\t")

EXON$Type<-c("exons")
INTRON$Type<-c("introns")
REP$Type<-c("repeats")


COMBO<-rbind(EXON,INTRON)
COMBO<-rbind(COMBO,REP)

eirplot<- ggplot(
  COMBO[COMBO$V1 != "MtDNA" , ],
  aes_string(
    x = "V2",
    y = COMBO[COMBO$V1 != "MtDNA" , ]$V7 * 100,
    color = "Type",
    fill = "Type",
    ordered = FALSE
  )
) + facet_grid(. ~ V1, scales  = "free") + theme_bw() + labs(title = "", x = "Genome Position (Mb)", y = "Fraction (%)") + theme(
  axis.text.x = element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.title = element_text(hjust = 0.5)
) + theme(legend.position = "top") + scale_colour_manual(
  values = c("#FF7C01", "#3bb3ad", "#696866"),
  name = ""
) + scale_fill_manual(
  values = c("#FF7C01", "#3bb3ad", "#696866"),
  name = ""
) + scale_x_continuous(
  labels = function(x)
    x / 1000000
) + theme(
  strip.background = element_rect(colour = "white", fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)  + geom_smooth(method = "loess",span = 0.4) + theme(strip.text = element_text(size = 10))




tiff(filename = "Fig1.tiff",width = 7.2,height = 3.2, res=300, units = "in")  ### Figure 2
plot(eirplot)
dev.off()


##############################################
################## Fig3 ######################
##############################################

library(ape)
#library(phytools)
library(ggtree)
library(patchwork)


data<-read.csv("C.brenneri_trees.txt",header=F,sep=" ", stringsAsFactors = F)
head(data)
data$CHR<-gsub("\\..*$","", perl=T,data$V1)
data$POS<-gsub(".*\\.([0-9]+)-.*$","\\1", perl=T,data$V1)


data_subset <- data[data$V2 != "" & !grepl(":[912]",data$V2) & !grepl(":0.[89]", data$V2) & !grepl("0.00000000",data$V2), ]

data_subset$heights <- sapply(data_subset$V2, function(tree_string) {
  tree <- read.tree(text = tree_string)
  max(node.depth.edgelength(tree))
})


#write.table(data_subset$V2, "ALL_trees2.nwk", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

TRSALL<-ape::read.tree("ALL_trees2.nwk")




set.seed(777)
SAMP<-sample(TRSALL,size=500)

sublist<-list("Cbren"="C.brenneri", "Csp48"="C.sp48", "Csp44"="C.sp44", "Csp51"="C.sp51", "Cdoug"="C.doughertyi", "Ctrop"="C.tropicalis", "Cwall"="C.wallacei", "Crema"="C.remanei", "Celeg"="C.elegans", "Csp54"="C.sp54")


SAMP2 <- lapply(SAMP, function(tree) {
  tree$tip.label <- Reduce(function(x, pattern) gsub(pattern, sublist[[pattern]], x), names(sublist), init = tree$tip.label)
  return(tree)
})
class(SAMP2)<- "multiPhylo"



treeALL <- ggtree(SAMP2,layout="slanted", color="#008080", alpha=.1)  + geom_tiplab(size=4, align=TRUE, linesize=.5,fontface = "italic")  + ggtitle("B") +
  theme(plot.title = element_text(face = "bold", hjust = 0))  + geom_treescale(x=0, y=0, width=0.1, color='black') + coord_cartesian(clip = "off") + theme(plot.margin = margin(0.2,2,0.2,0.2, "cm"))



treeALL




bra<-read.csv("C.brenneri_branches.txt",header=F, sep=" ")
head(bra)
levels(as.factor(bra$V1))
table((bra$V1))
bra<-bra[bra$V1 %in% c("I", "II", "III", "IV", "V", "X"),]
bra$V1<-as.character(bra$V1)
bra<-data.frame(bra)
bra$Type <- "C.brenneri&C.sp48\nbranch"
bra<-bra[,-3]
colnames(bra)<-c("CHR",     "POS",     "heights", "Type")


hei<-data_subset[data_subset$CHR!="MtDNA",c(3,4,5)]
hei$Type <- "Total tree height"
hei$POS<-as.numeric(as.character(hei$POS))



phastcons <- read.csv("Cbren_combo_phastcons_50kb.wig",header=F, sep=" ")
phastcons$CHR<-gsub("\\..*$","", perl=T,phastcons$V1)
phastcons$POS<-gsub(".*\\.([0-9]+)-.*$","\\1", perl=T,phastcons$V1)
phastcons$Type <- "phastCons score"
phastcons <-phastcons[,c(3,4,2,5)]
colnames(phastcons)<-c("CHR",     "POS",     "heights", "Type")
phastcons$POS<-as.numeric(as.character(phastcons$POS))
phastcons$heights <- as.numeric(as.character(phastcons$heights))

hei<-rbind(hei,bra)
hei<-rbind(hei,phastcons)


hei$Type <- factor(hei$Type, levels = c("Total tree height","C.brenneri&C.sp48\nbranch","phastCons score"))

hei <- hei[complete.cases(hei),]

#table(hei[hei$heights > 1 & !is.na(hei$heights),]$Type)

#Total tree height Cbren&C.sp48\nbranch      phastCons score
#16                   15                    0
heights<-ggplot(hei[hei$heights < 1,], aes(x = POS, y = heights, ordered = FALSE))  + facet_grid(Type ~ CHR, scales =
                                                                                   "free") + theme_bw(base_size = 12) + labs(title = "C", x = "Genome Position (Mb)", y = "value") + theme(
                                                                                     axis.text.x = element_text(
                                                                                       angle = 0,
                                                                                       vjust = 0.5,
                                                                                       hjust = 0.5
                                                                                     ),
                                                                                     plot.title = element_text(face = "bold", hjust = 0)
                                                                                   ) + theme(legend.position = "none") + theme(
                                                                                     strip.background = element_blank(),
                                                                                     panel.grid.major = element_blank(),
                                                                                     panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))
                                                                                   )  + geom_point(size = 0.45, alpha = 0.35) + geom_smooth(
                                                                                     method = "loess",
                                                                                     span = 0.4,
                                                                                     color = "#3bb3ad",
                                                                                     linetype = "solid",
                                                                                     se = F,
                                                                                     size = 1.1
                                                                                   ) + theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12)) + scale_x_continuous(
                                                                                     labels = function(x)
                                                                                       x / 1000000
                                                                                   ) + coord_cartesian(clip = "off")













### topology

topl<-read.csv("C.brenneri_topol.txt",header=F, sep=" ")
topl$Type<-"Other topology"
topl[topl$V4=="",]$Type <- "NA"
topl$Type[grepl("(Cbren,Csp48)",topl$V4)] <- "((C.brenneri,C.sp48),others);"

topl <-topl[topl$V1!="MtDNA",]
topl<-as.data.frame(topl)
topl$Type <-factor(topl$Type,levels=c("((C.brenneri,C.sp48),others);","Other topology","NA"))


#library(data.table)
##ff7f00
toplot<-  ggplot(topl, aes(y = V2))  +
  coord_flip() +   facet_grid(.~V1,scales="free") + geom_rect(data=topl,aes(ymin=V2, ymax=V3, xmin=10, xmax=20,fill=Type,color=Type)) + theme_bw() +  theme(
    strip.background = element_rect(colour = "white", fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),legend.title=element_blank()
  )  + scale_y_continuous(
    labels = function(x)
      x / 1000000
  ) + scale_colour_manual(values = c("#3f918a","#ff7f00","grey95"), name="",
                          labels = c(
                            parse(text = expression(paste("((", italic(C.brenneri), ",", italic(C.sp48), "),others);"))),
                            parse(text = "plain('Other topology')"),
                            parse(text = "plain('NA')")
                          )) + scale_fill_manual(values = c("#3f918a","#ff7f00","grey95"),name="",
                                                 labels = c(
                                                   parse(text = expression(paste("((", italic(C.brenneri), ",", italic(C.sp48), "),others);"))),
                                                   parse(text = "plain('Other topology')"),
                                                   parse(text = "plain('NA')")
                                                 )) + labs(title = "A", y = "Genome Position (Mb)", x = "Topology")  + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(face = "bold", hjust = 0), legend.text = element_text(hjust = 0)) + theme(legend.position = "top")






FIG2<- (toplot/treeALL/heights) +
  plot_layout(heights  = c(1,8,15))


tiff(filename = "Fig2.tiff",width = 7.2,height = 10.9, res=300, units = "in") ###Figure 3
plot(FIG2)
dev.off()


##############################################
################## Fig4 ######################
##############################################

##for figure 4 see file comparative_genomics/orthologs/4_run_genespace_get_Fig4.R



##############################################
################## Fig5 ######################
##############################################

library(patchwork)
library(reshape2)

#a species tree from Orthofinder2
TREESP<-c("(Cinop:0.091484,(Celeg:0.134144,((Crema:0.127918,(Cnigo:0.0504089,Cbrig:0.0334806)0.90561:0.117991)0.440453:0.026276,(Ctrop:0.129409,Cbren:0.125876)0.533563:0.0320385)0.636122:0.036656)1:0.091484);")

TREESP <- read.tree(text = TREESP)
dist_matrix <- cophenetic(TREESP)



#a tree in a nice order with full names
TREESP2<-c("(C. inopinata:0.091484,(C. elegans:0.134144,((C. remanei:0.127918,(C. nigoni:0.0504089,C. briggsae:0.0334806)0.90561:0.117991)0.440453:0.026276,(C. tropicalis:0.129409,C. brenneri:0.125876)0.533563:0.0320385)0.636122:0.036656)1:0.091484);")
TREESP2 <- read.tree(text = TREESP2)
TREESPROT<-rotateConstr(TREESP2,c("C.inopinata","C.elegans","C.tropicalis","C.brenneri","C.briggsae","C.nigoni","C.remanei"))

plot(TREESPROT)

OGcounts<-read.csv("Orthogroups.GeneCount.tsv",sep='\t')
OGcounts_filt<-OGcounts[,c(2:8)]
#remove singletons
OGcounts_filt<-OGcounts_filt[which(rowSums(OGcounts_filt > 0) > 1),]
corOG<-cor(OGcounts_filt)



# Melt the lower triangle of the first matrix
melted1 <-melt(dist_matrix)
#dist_matrix2<- dist_matrix
#dist_matrix2[upper.tri(dist_matrix2,diag = TRUE)] <- NA
corOG2<-corOG
corOG2[upper.tri(corOG2,diag = TRUE)] <- NA
#melted1.2 <- melt(dist_matrix2,na.rm = TRUE)
#melted2 <-melt(corOG2,na.rm = TRUE)
lower<-melt(corOG2,na.rm = TRUE)[,c(1,2)]
lower$sp1 <- "Selfing"
lower[lower$Var1 %in% c("Cbren","Cinop","Crema", "Cnigo"),]$sp1 <- "Outcrossing"
lower$sp2 <- "Selfing"
lower[lower$Var2 %in% c("Cbren","Cinop","Crema", "Cnigo"),]$sp2 <- "Outcrossing"
lower$comp <- paste0(lower$sp1,"/",lower$sp2)
swap_condition <- lower$comp == "Outcrossing/Selfing"
lower[swap_condition, c("Var1", "Var2")] <- lower[swap_condition, c("Var2", "Var1")]
lower <- lower[,c(1,2)]




melted2 <-melt(corOG) # let's use r not r2
#melted1 <- subset(melted1, Var1 != Var2)
#melted2 <- subset(melted2, Var1 != Var2)

OGCORandDIST2 <- merge(melted1, melted2, by = c("Var1", "Var2"), all.y = TRUE)
OGCORandDIST2 <- merge(OGCORandDIST2,lower,by = c("Var1", "Var2"), all.y = TRUE)

OGCORandDIST2$sp1 <- "Selfing"
OGCORandDIST2[OGCORandDIST2$Var1 %in% c("Cbren","Cinop","Crema", "Cnigo"),]$sp1 <- "Outcrossing"
OGCORandDIST2$sp2 <- "Selfing"
OGCORandDIST2[OGCORandDIST2$Var2 %in% c("Cbren","Cinop","Crema", "Cnigo"),]$sp2 <- "Outcrossing"

OGCORandDIST2$comp <- paste0(OGCORandDIST2$sp1,"/",OGCORandDIST2$sp2)




pltre<-ggtree(TREESPROT,layout="slanted", color="#3f918a", alpha=0.75,size=1.5,ladderize=F) + geom_tiplab(fontface="italic",size=3.5, align=FALSE, linesize=.5)  +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  + geom_treescale(x=0, y=0, width=0.1, fontsize=3.5, color='grey40') + theme(plot.margin = margin(0.2,0.7,0.2,0.2, "cm")) + ggplot2::xlim(0, 0.4) +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data)
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = FALSE, params = list(grob = grob,
                                           xmin = xmin, xmax = xmax,
                                           ymin = ymin, ymax = ymax))
}


OGCORandDIST2$comp <- as.character(OGCORandDIST2$comp)
corOGpl<- ggplot(OGCORandDIST2, aes(x = value.x, y = value.y, ordered = FALSE, shape = comp)) + theme_bw(base_size = 12) + labs(
  title = "A",
  x = "Genetic distance",
  y = "Correlation in OG sizes (r)") + theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    plot.title = element_text(hjust = 0.5)
  ) + theme(legend.position = "right") + scale_colour_manual(values = c("#A49C74")) + scale_fill_manual(values = c("#A49C74"))  + theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = rel(1))
  )  + geom_point(size = 3, alpha = 1) + theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12)) + annotation_custom2(grob=ggplotGrob(pltre), data=TREESPROT,
                                                                                                                                                         xmin = 0.18, xmax=0.465, ymin=0.27, ymax=0.49)   +
  theme(plot.title = element_text(face = "bold", hjust = 0)) +
  theme(plot.margin = margin(0.2,0.7,0.2,0.5, "cm")) +
  geom_label(aes(x = 0.13, y = 0.07, label = "18,282 OGs"), fill = "#ffe0a6",label.size = NA) +
  scale_shape_manual(values = c("Outcrossing/Outcrossing" = 15, "Selfing/Outcrossing" = 7, "Selfing/Selfing" = 4),name="")










corOGpl


#a new file with info 0/1, sizes, phase, splicing sites, 7238 OG, all genes have at leat one intron (90 OG where removed)
intstrALL<- read.csv("7sp_all_intron_info_ALL.FIN.tab",header=F,sep="\t")
#"Cbren" "Cbrig" "Celeg" "Cinop" "Cnigo" "Crema" "Ctrop"

prefixes <- c("Cbren", "Cbrig", "Celeg", "Cinop", "Cnigo", "Crema", "Ctrop")
suffixes <- c("intr", "size", "phase", "ss")  # Add other suffixes as needed

# Create vector of new column names
new_column_names <- c(paste0(rep(prefixes, each = length(suffixes)), "_", rep(suffixes, times = length(prefixes))),"Orthogroup")

colnames(intstrALL) <- new_column_names
#intstrALL<-intstrALL[complete.cases(intstrALL),] #all complete:)

intrBIN <- as.matrix(intstrALL[, grep("_intr", names(intstrALL))])

#remove singletons
intrBIN2<-intrBIN[which(rowSums(intrBIN > 0) > 1),]
corBIN<-cor(intrBIN2)
#cor(intrBIN)

melted4 <-melt(corBIN) #r not r^2
melted4$Var1 <- gsub("_intr","",melted4$Var1)
melted4$Var2 <- gsub("_intr","",melted4$Var2)
melted4 <- subset(melted4, Var1 != Var2)




INTRONBINandDIST <- merge(melted1, melted4, by = c("Var1", "Var2"), all.x = TRUE)

INTRONBINandDIST2 <- merge(INTRONBINandDIST,lower,by = c("Var1", "Var2"), all.y = TRUE)
INTRONBINandDIST2$sp1 <- "Selfing"
INTRONBINandDIST2[INTRONBINandDIST2$Var1 %in% c("Cbren","Cinop","Crema", "Cnigo"),]$sp1 <- "Outcrossing"
INTRONBINandDIST2$sp2 <- "Selfing"
INTRONBINandDIST2[INTRONBINandDIST2$Var2 %in% c("Cbren","Cinop","Crema", "Cnigo"),]$sp2 <- "Outcrossing"

INTRONBINandDIST2$comp <- paste0(INTRONBINandDIST2$sp1,"/",INTRONBINandDIST2$sp2)


rrr<- rowSums(intrBIN)

intron_spec <- as.data.frame(table(rrr))

barpl <-ggplot(intron_spec, aes(x = rrr, y = Freq, ordered = FALSE)) + theme_bw(base_size = 10) + labs(
  title = "",
  x = "Number of species",
  y = "Common introns") + theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    plot.title = element_text(hjust = 0.5)
  ) + theme(legend.position = "none") + theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = rel(1))
  )  + geom_bar(stat="identity", col="#3f918a", fill="#3f918a") + theme(strip.text = element_text(size = 10))  + theme(strip.text = element_text(size = 10), panel.border = element_blank(), axis.line = element_line(colour = "grey40")) +  scale_y_continuous(sec.axis = sec_axis(~ . /sum(intron_spec$Freq)*100, name = "%"))







corINTpl2<- ggplot(INTRONBINandDIST2, aes(x = value.x, y = value.y, ordered = FALSE,shape=comp)) + theme_bw(base_size = 12) + labs(
  title = "B",
  x = "Genetic distance",
  y = "Correlation in gene structures (r)") + theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    plot.title = element_text(hjust = 0.5)
  ) + theme(legend.position = "none") + scale_colour_manual(values = c("#A49C74")) + scale_fill_manual(values = c("#A49C74"))  + theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = rel(1))
  )  + geom_point(size = 3, alpha = 1) + theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12))   + theme(plot.title = element_text(face = "bold", hjust = 0)) + theme(plot.margin = margin(0.2,0.7,0.2,0.5, "cm")) + geom_label(aes(x = 0.13, y = -0.07, label = "7,238 1-to-1 OGs\n167,873 introns*"), fill = "#ffe0a6",alpha=0.5,label.size = NA) + annotation_custom2(grob=ggplotGrob(barpl), data=intron_spec,
                                                                                                                                                                                                                                                                                                                                                                                                             xmin = 0.2, xmax=0.435, ymin=0.25, ymax=0.65) + geom_hline(yintercept = 0, linetype = "dotted", linewidth=0.9, color = "grey80") +
  scale_shape_manual(values = c("Outcrossing/Outcrossing" = 15, "Selfing/Outcrossing" = 7, "Selfing/Selfing" = 4),name="")








corINTpl2




FIG4<- (corOGpl +  corINTpl2 + plot_layout(guides = 'collect')  & theme(legend.position = 'bottom'))



tiff(filename = "Fig4.2.tiff",width = 12.2,height = 6, res=300, units = "in") # Figure 5
plot(FIG4)
dev.off()



##############################################
################## Fig6 ######################
##############################################



my_filesS<-list.files(pattern="_gt_stat_")
my_names<-gsub(".txt","", perl=T, my_filesS)

STAT<-c();for (i in 1:(length(my_filesS)) ) { B=read.csv(file=my_filesS[i], header=F, sep=" "); STAT<-rbind(STAT,data.frame(B,sample=my_names[i]));}

STAT$V1<-gsub(":","",STAT$V1)

STAT$Species <- gsub("([a-zA-Z]+)(_.*)","\\1",perl=T,STAT$sample)
STAT$Type <- gsub("(.*)_([a-zA-Z]+)$","\\2",perl=T,STAT$sample)

STAT$Species <-gsub("Cbren","C. brenneri", STAT$Species)
STAT$Species <-gsub("Celeg","C. elegans", STAT$Species)
STAT$Species <-gsub("Crema","C. remanei", STAT$Species)
STAT$Species <-gsub("Cbrig","C. briggsae",STAT$Species)
STAT$Species <-gsub("Cinop","C. inopinata",STAT$Species)
STAT$Species <-gsub("Ctrop", "C. tropicalis",STAT$Species)
STAT$Species <-gsub("Cnigo", "C. nigoni",STAT$Species)
#"Cbren" "Cbrig" "Celeg" "Cinop" "Cnigo" "Crema" "Ctrop"

STAT$Species <- factor(STAT$Species, levels = c("C. brenneri","C. remanei","C. nigoni","C. inopinata","C. briggsae","C. tropicalis","C. elegans"))
#"exonlength"   "exonnumber"   "genelength"   "intronlength"

STAT$Type <- gsub("exonlength", "exon length", STAT$Type)
STAT$Type <- gsub("exonnumber", "exon number", STAT$Type)
STAT$Type <- gsub("intronlength", "intron length", STAT$Type)
STAT$Type <- gsub("genelength", "gene length", STAT$Type)
STAT$Type <- factor(STAT$Type, levels = c("gene length", "exon length", "intron length", "exon number"))

STAT$Gene<-gsub("[a-zA-Z]+_([a-zA-Z-0-9._]+)_gt_stat_.*$","\\1",perl=T,STAT$sample)
STAT$Gene <- gsub("justall", "all", STAT$Gene)

#STAT2$Gene<-as.character(STAT2$Gene)
STAT$Gene <- factor(STAT$Gene, levels = c("1-to-1", "all"))

STAT<-as.data.frame(STAT)
STAT$V1<-as.numeric(as.character(STAT$V1))
STAT$V2<- as.numeric(as.character(STAT$V2))
STAT<-STAT[complete.cases(STAT),]
STAT3 <- data.frame(value = rep(STAT$V1, STAT$V2),
                    Species = rep(STAT$Species, STAT$V2),
                    Type = rep(STAT$Type, STAT$V2),
                    Gene = rep(STAT$Gene, STAT$V2))



STAT3med <- STAT3[STAT3$Type %in% c("exon length", "intron length") &
                    STAT3$Gene %in% c("1-to-1", "all"),] %>% group_by(Species, Gene,Type) %>%
  summarize(med = median(value, na.rm = TRUE), mean=mean(value, na.rm = TRUE))



exinlen <-  ggplot(STAT3[STAT3$Type %in% c("exon length", "intron length") &
                           STAT3$Gene %in% c("1-to-1", "all"),],
                   aes(
                     y = value,
                     x = Species,
                     fill = Species,
                     color = Species,
                     ordered = FALSE
                   )) +
  facet_grid(Gene ~ Type, scales = "free") +
  scale_color_manual(values =
                       c("black", "#737170", "#3bb3ad", "#4fa0b8", "#FFA500", "#ff7f00", "#f0b326")) +
  scale_fill_manual(values =
                      c("black", "#737170", "#3bb3ad", "#4fa0b8", "#FFA500", "#ff7f00", "#f0b326")) +
  geom_violin(position = "identity", alpha = 0.4) + ylim(0, 500) + coord_cartesian(clip ="off") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    plot.title = element_text(hjust = 0, face = "bold"),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = rel(1))
  ) +
  geom_point(data = STAT3med, aes(x = Species, y = med, color = Species), shape = 23, size = 2, color = "black", fill = "white") +
  geom_point(data = STAT3med, aes(x = Species, y = mean, color = Species), size = 2,shape = 4,color = "black", fill = "white") +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.position = "none") +
  coord_flip(clip = "off") +
  labs(y = "bp", title = "A")




### as a histogram
STAT4<- STAT3[STAT3$Type=="exon number" & STAT3$Gene %in% c("1-to-1", "all"),]

STAT4group <- STAT4 %>%
  group_by(Species, Gene, value) %>%
  summarise(Fraction = n() / nrow(STAT4))

#normalize
STAT4group <- STAT4group %>%
  group_by(Species, Gene) %>%
  mutate(Fraction = Fraction / sum(Fraction))

STAT4group<- as.data.frame(STAT4group)
STAT4med <- STAT4 %>% group_by(Species, Gene,Type) %>%
  summarize(med = median(value, na.rm = TRUE), mean=mean(value, na.rm = TRUE))
STAT4group$Species <- factor(STAT4group$Species, levels = rev(levels(STAT4group$Species)))


numpl3<- ggplot(STAT4group,aes(x=value, y=Fraction*100, fill = Species, color = Species, ordered = FALSE))  +    facet_grid(Species ~ Gene)   +scale_color_manual(values=rev(c("black","#737170","#3bb3ad","#4fa0b8", "#FFA500","#ff7f00","#f0b326")),name="") + scale_fill_manual(values=rev(c("black","#737170","#3bb3ad","#4fa0b8", "#FFA500","#ff7f00","#f0b326")),name="") +
  geom_col(position = "identity", alpha=0.4) + theme_bw(base_size = 12) +coord_flip(clip="off") +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) + coord_cartesian(clip="off") + xlim(1,16)  + labs(title="B", x = "Number of exons", y="% of genes") + theme(legend.position = "none") + theme(legend.text = element_text(face = "bold.italic")) +
  geom_vline(data = STAT4med, aes(xintercept = med), linetype = "dotted") +
  geom_vline(data = STAT4med, aes(xintercept = mean), linetype = "dashed") +
  theme(strip.text.y = element_text(face = "italic",angle = 0,hjust=0))



FIG5 <- (exinlen/numpl3) + plot_layout(heights  = c(1,1.5))

tiff(filename = "Fig5.tiff",width = 6.8,height = 8.7, res=300, units = "in")  ### Figure 6
plot(FIG5)
dev.off()


##############################################
################## FigS7,FigS8 ###############
##############################################

###A. gene lengths

STAT3medGENE <- STAT3[STAT3$Type %in% c("gene length") &
                        STAT3$Gene %in% c("1-to-1", "all"),] %>% group_by(Species, Gene,Type) %>%
  summarize(med = median(value, na.rm = TRUE), mean=mean(value, na.rm = TRUE))




genelen <-  ggplot(STAT3[STAT3$Type %in% c("gene length") &
                           STAT3$Gene %in% c("1-to-1", "all"),],
                   aes(
                     y = value,
                     x = Species,
                     fill = Species,
                     color = Species,
                     ordered = FALSE
                   )) +
  facet_grid(Gene ~ Type, scales = "free") +
  scale_color_manual(values =
                       c("black", "#737170", "#3bb3ad", "#4fa0b8", "#FFA500", "#ff7f00", "#f0b326")) +
  scale_fill_manual(values =
                      c("black", "#737170", "#3bb3ad", "#4fa0b8", "#FFA500", "#ff7f00", "#f0b326")) +
  geom_violin(position = "identity", alpha = 0.4)  + coord_cartesian(clip ="off") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    plot.title = element_text(hjust = 0, face = "bold"),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = rel(1))
  ) +
  geom_point(data = STAT3medGENE, aes(x = Species, y = med, color = Species), shape = 23, size = 2, color = "black", fill = "white") +
  geom_point(data = STAT3medGENE, aes(x = Species, y = mean, color = Species), size = 2,shape = 4,color = "black", fill = "white") +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.position = "none") +
  coord_flip(clip = "off") +
  labs(y = "bp", title = "") + ylim(0, 5000)





## number of introns

library(stats)
library(dplyr)
library(esvis)
library(tidyr)
library(data.table)
library(reshape2)
library(corrplot)
library(RColorBrewer)
library(ggplotify)
library(gridGraphics)

sp_real_order <- c("C. brenneri","C. briggsae","C. elegans","C. inopinata","C. nigoni","C. remanei","C. tropicalis")
sp_phylo_order <- c("C. inopinata","C. elegans","C. tropicalis","C. brenneri","C. briggsae","C. nigoni","C. remanei")

#change margins
par(mar = c(1, 1, 1, 1))


pairwise_resultEX <- pairwise.wilcox.test(STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="exon number",]$value, STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="exon number",]$Species,na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE, p.adjust.method = "bonferroni")



grouped_dataEX <- STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="exon number",] %>%
  group_by(Species)

grouped_dataEX %>% coh_d(value~ Species) -> chdmaEX


wide_dataEX <- dcast(chdmaEX, Species_ref ~ Species_foc, value.var = "coh_d")
wide_dataEX <- wide_dataEX[, -1]
rownames(wide_dataEX) <-colnames(wide_dataEX)
wide_dataEX<-as.matrix(wide_dataEX)

symmetric_matrixpvalEX <- pairwise_resultEX$p.value
first_rowEX <- setNames(rep(NA, ncol(symmetric_matrixpvalEX)), colnames(symmetric_matrixpvalEX))
final_matrixEX <- rbind(first_rowEX, symmetric_matrixpvalEX)
last_colEX <- setNames(rep(NA, nrow(final_matrixEX)), rownames(final_matrixEX))
final_matrixEX <-cbind(final_matrixEX,last_colEX)
names2EX<-c(colnames(final_matrixEX)[1:6],"C. elegans")
colnames(final_matrixEX) <-names2EX
rownames(final_matrixEX) <-names2EX
final_matrixEX[upper.tri(final_matrixEX)] <- t(final_matrixEX)[upper.tri(final_matrixEX)]
final_matrixEX<-as.matrix(as.data.frame(final_matrixEX))

corrplot(wide_dataEX[sp_phylo_order,sp_phylo_order], insig = "label_sig",col.lim = c(-0.2, 0.2),sig.level = c(.0001, .001, .01), mar=c(0,0,1,0), title="B. Exon number                                                                                                ",  is.corr = FALSE,p.mat=final_matrixEX[sp_phylo_order, sp_phylo_order],pch.cex = 0.8, diag = FALSE, font = 3, tl.col = 'black', col=colorRampPalette(c("#633f01","white","#016163"))(8),cl.ratio = 0.25, type = "upper")


cor_plotEX <- grab_grob()
cor_plotEX <- as.ggplot(cor_plotEX)

cor_plotEX  <- cor_plotEX  +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))




####lenght of introns

pairwise_result <- pairwise.wilcox.test(STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="intron length",]$value, STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="intron length",]$Species, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE, p.adjust.method = "bonferroni")

grouped_data <- STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="intron length",] %>%
  group_by(Species)


grouped_data %>% coh_d(value~ Species) -> chdma


wide_data <- dcast(chdma, Species_ref ~ Species_foc, value.var = "coh_d")
wide_data <- wide_data[, -1]
rownames(wide_data) <-colnames(wide_data)
wide_data<-as.matrix(wide_data)

symmetric_matrixpval <- pairwise_result$p.value
first_row <- setNames(rep(NA, ncol(symmetric_matrixpval)), colnames(symmetric_matrixpval))
final_matrix <- rbind(first_row, symmetric_matrixpval)
last_col <- setNames(rep(NA, nrow(final_matrix)), rownames(final_matrix))
final_matrix <-cbind(final_matrix,last_col)
names2<-c(colnames(final_matrix)[1:6],"C. elegans")
colnames(final_matrix) <-names2
rownames(final_matrix) <-names2
final_matrix[upper.tri(final_matrix)] <- t(final_matrix)[upper.tri(final_matrix)]


corrplot(wide_data[sp_phylo_order,sp_phylo_order], col.lim = c(-0.3, 0.3), insig = "label_sig",sig.level = c(.0001, .001, .01), mar=c(0,0,1,0), title="A. Intron lenght                                                                                                ",  is.corr = FALSE,p.mat=final_matrix[sp_phylo_order, sp_phylo_order],pch.cex = 0.8, diag = FALSE, font = 3, tl.col = 'black', col=colorRampPalette(c("#633f01","white","#016163"))(6),cl.ratio = 0.25, type = "upper")


cor_plotLEN <- grab_grob()
cor_plotLEN <- as.ggplot(cor_plotLEN)
cor_plotLEN <- cor_plotLEN +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))



#### length of genes
pairwise_resultG <- pairwise.wilcox.test(STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="gene length",]$value, STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="gene length",]$Species, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE, p.adjust.method = "bonferroni")

grouped_dataG <- STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="gene length",] %>%
  group_by(Species)


grouped_dataG %>% coh_d(value~ Species) -> chdmaG


wide_dataG <- dcast(chdmaG, Species_ref ~ Species_foc, value.var = "coh_d")
wide_dataG <- wide_dataG[, -1]
rownames(wide_dataG) <-colnames(wide_dataG)
wide_dataG<-as.matrix(wide_dataG)

symmetric_matrixpvalG <- pairwise_resultG$p.value
first_rowG <- setNames(rep(NA, ncol(symmetric_matrixpvalG)), colnames(symmetric_matrixpvalG))
final_matrixG <- rbind(first_rowG, symmetric_matrixpvalG)
last_colG <- setNames(rep(NA, nrow(final_matrixG)), rownames(final_matrixG))
final_matrixG <-cbind(final_matrixG,last_colG)
names2G<-c(colnames(final_matrixG)[1:6],"C. elegans")
colnames(final_matrixG) <-names2G
rownames(final_matrixG) <-names2G
final_matrixG[upper.tri(final_matrixG)] <- t(final_matrixG)[upper.tri(final_matrixG)]


corrplot(wide_dataG[sp_phylo_order,sp_phylo_order], insig = "label_sig",col.lim = c(-0.4, 0.4),sig.level = c(.0001, .001, .01), mar=c(0,0,1,0), title="C. Gene lenght                                                                                                ",  is.corr = FALSE,p.mat=final_matrixG[sp_phylo_order, sp_phylo_order],pch.cex = 0.8, diag = FALSE, font = 3, tl.col = 'black', col=colorRampPalette(c("#633f01","white","#016163"))(8),cl.ratio = 0.25, type = "upper")


cor_plotLENG <- grab_grob()
cor_plotLENG <- as.ggplot(cor_plotLENG)
cor_plotLENG <- cor_plotLENG +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))




FIGS3 <- (cor_plotLEN + cor_plotEX + cor_plotLENG)

tiff(filename = "FigS3B.tiff",width = 22.6,height = 9.6, res=300, units = "in") ## Figure S8

FIGS3
dev.off()

tiff(filename = "FigS3A.tiff",width = 6.8,height = 3.48, res=300, units = "in") ## Figure S7
genelen
dev.off()



#############################################


##############################################
################## FigS5, FigS9 ##############
##############################################

library(rcompanion)
library(corrplot)
library(ggplotify)
library(gridGraphics)


intrones_SIZE <- as.matrix(intstrALL[, grep("_size", names(intstrALL))])
cor_matrixSIZE <- matrix(NA, nrow = ncol(intrones_SIZE), ncol = ncol(intrones_SIZE))
cor_matrixSIZECOR <- matrix(NA, nrow = ncol(intrones_SIZE), ncol = ncol(intrones_SIZE))
cor_matrixSIZETEST <- matrix(NA, nrow = ncol(intrones_SIZE), ncol = ncol(intrones_SIZE))




calculate_nonzero_mean_and_wilx <- function(col1, col2) {
  nonzero_rows <- which(col1 > 0 & col2 > 0)
  data <- as.data.frame(cbind(col1[nonzero_rows],col2[nonzero_rows]))
  data <- data[(data$V1 > quantile(data$V1, probs=c(.05, .95), na.rm = FALSE)[1] & data$V1 < quantile(data$V1, probs=c(.05, .95), na.rm = FALSE)[2] ) &
                 (data$V2 > quantile(data$V2, probs=c(.05, .95), na.rm = FALSE)[1] & data$V2 < quantile(data$V2, probs=c(.05, .95), na.rm = FALSE)[2] ),]
  if (length(nonzero_rows) > 0) {
    mn<-mean(data$V1-data$V2)
    wt <- wilcox.test(data$V1, data$V2, paired = TRUE)$p.value
    cor <- cor(data$V1, data$V2)

    result <- list(mean = mn, pval = wt,cor=cor)
    return(result)

  } else {
    result <- list(mean = NA, pval = NA,cor=NA)
    return(result)
  }
}



calculate_nonzero_mean_and_wilx2 <- function(col1, col2) {
  nonzero_rows <- which(col1 > 0 & col2 > 0)
  data <- as.data.frame(cbind(col1[nonzero_rows],col2[nonzero_rows]))
  #data <- data[(data$V1 > quantile(data$V1, probs=c(.05, .95), na.rm = FALSE)[1] & data$V1 < quantile(data$V1, probs=c(.05, .95), na.rm = FALSE)[2] ) &
  #               (data$V2 > quantile(data$V2, probs=c(.05, .95), na.rm = FALSE)[1] & data$V2 < quantile(data$V2, probs=c(.05, .95), na.rm = FALSE)[2] ),]
  if (length(nonzero_rows) > 0) {
    #num <- nrow(data)
    #perc <-
    mn<-mean(data$V1-data$V2)
    rat <-mean(data$V1/data$V2)
    wt <- wilcox.test(data$V1, data$V2, paired = TRUE)$p.value
    cor <- cor(data$V1, data$V2)

    result <- list(mean = mn, pval = wt,cor=cor,num=num)
    return(result)

  } else {
    result <- list(mean = NA, pval = NA,cor=NA)
    return(result)
  }
}


calculate_nonzero_num_and_percent <- function(col1, col2) {
  col1<- as.numeric(as.character(col1))
  col2<- as.numeric(as.character(col2))
  nonzero_rows <- which(col1 > 0 & col2 > 0)

  if (length(nonzero_rows) > 0) {
    num <- length(nonzero_rows)
    perc1 <- length(nonzero_rows)/length(which(col1 > 0))*100
    perc2 <- length(nonzero_rows)/length(which(col2 > 0))*100
    tot1 <- sum(col1[nonzero_rows])
    tot2 <- sum(col2[nonzero_rows])

    result <- list(num=num,perc1=perc1,perc2=perc2,tot1=tot1,tot2=tot2 )
    return(result)

  } else {
    result <- list(num=NA,perc1=NA,perc2=NA,tot1=NA,tot2=NA)
    return(result)
  }
}


calculate_nonzero_num_and_percent_diag <- function(col1, df) {
  col1<- as.numeric(as.character(col1))
  #col2<- as.numeric(as.character(col2))
  #OGcounts_filt[which(rowSums(df > 0) > 1),]

  nonzero_rows <- which(col1 > 0 & rowSums(df > 0) == 1)

  if (length(nonzero_rows) > 0) {
    perc1 <- length(nonzero_rows)/length(col1)*100
    tot1 <- sum(col1[nonzero_rows])

    result <- list(perc1=perc1,tot1=tot1 )
    return(result)

  }

}




cor_matrixOGPER <- matrix(NA, nrow = ncol(OGcounts_filt), ncol = ncol(OGcounts_filt))
cor_matrixOGSUM <- matrix(NA, nrow = ncol(OGcounts_filt), ncol = ncol(OGcounts_filt))

OGcounts2<- OGcounts[,c(2:8)]

for (i in 1:(ncol(OGcounts2))) {
  for (j in (i):ncol(OGcounts2)) {
    if (i != j) {
      res <- calculate_nonzero_num_and_percent(OGcounts2[, i], OGcounts2[, j])
      cor_matrixOGPER[i, j] <- res$perc1
      cor_matrixOGPER[j, i] <- res$perc2
      cor_matrixOGSUM[i, j] <- res$tot1
      cor_matrixOGSUM[j, i] <- res$tot2
    } else {
      res <- calculate_nonzero_num_and_percent_diag(OGcounts2[, i], OGcounts2)
      cor_matrixOGPER[i, j] <- res$perc1
      cor_matrixOGSUM[i, j] <- res$tot1
    }
  }
}



colnames(cor_matrixOGPER) <- sp_real_order
colnames(cor_matrixOGSUM) <- sp_real_order
rownames(cor_matrixOGPER) <- sp_real_order
rownames(cor_matrixOGSUM) <- sp_real_order

cor_matrixOGPER2 <- cor_matrixOGPER[sp_phylo_order, sp_phylo_order]
cor_matrixOGSUM2<- cor_matrixOGSUM[sp_phylo_order, sp_phylo_order]


cor_matrixOGPER2
#C. inopinata C. elegans C. tropicalis C. brenneri C. briggsae C. nigoni C. remanei
#C. inopinata      1.472563  93.338859     90.950192   88.192299  93.7608319 93.873860  89.254766
#C. elegans       86.634494   1.403103     90.047559   88.012309  94.4887397 92.761225  89.355155
#C. tropicalis    81.647839  87.093283      1.861542   87.810323  89.9952648 91.179057  88.094433
#C. brenneri      77.249026  83.057224     85.677513    4.088909  86.5289420 87.472774  84.991090
#C. briggsae      75.687348  82.177616     80.924574   79.744526   0.3797175 96.155718  81.697080
#C. nigoni        72.015723  76.669172     77.917799   76.611365  91.3810047  2.486687  79.599977
#C. remanei       76.513145  82.526969     84.122473   83.179381  86.7579614 88.947742   3.648993

cor_matrixOGSUM2 #on the diagonal # of OGs unique for the speceis
#C. inopinata C. elegans C. tropicalis C. brenneri C. briggsae C. nigoni C. remanei
#C. inopinata          3223      14714         14235       13475       14506     15279      13799
#C. elegans           14776       1483         15838       15540       16505     16230      15912
#C. tropicalis        14307      15348          2157       16018       15967     16227      16076
#C. brenneri          14212      15592         17137        4480       16685     17152      17227
#C. briggsae          14360      15850         15909       15757         246     19260      16198
#C. nigoni            15383      16638         17797       17499       22872      2367      18499
#C. remanei           14563      16176         17148       17217       17409     18220       3830

corrplot(cor_matrixOGPER2,  mar=c(0,0,3,0), col.lim = c(0, 100),title="B. % of shared OGs (species-specific on the diagonal)                                                                        ",  is.corr = FALSE,  pch.cex = 0.8,font = 3, tl.col = 'black',cl.ratio = 0.25,col=colorRampPalette(c("#633f01","white","#016163"))(10), p.mat = cor_matrixOGSUM2, insig = "p-value", sig.level = -1)

grab_grob <- function(){
  grid.echo()
  grid.grab()
}


cor_plotOGPER <- grab_grob()
cor_plotOGPER <- as.ggplot(cor_plotOGPER)
cor_plotOGPER  <- cor_plotOGPER  +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))


#### OG sizes
cor_matrixOGSIZE <- matrix(NA, nrow = ncol(OGcounts_filt), ncol = ncol(OGcounts_filt))
cor_matrixOGSIZECOR <- matrix(NA, nrow = ncol(OGcounts_filt), ncol = ncol(OGcounts_filt))
cor_matrixOGSIZETEST <- matrix(NA, nrow = ncol(OGcounts_filt), ncol = ncol(OGcounts_filt))

for (i in 1:(ncol(OGcounts_filt) - 1)) {
  for (j in (i + 1):ncol(OGcounts_filt)) {
    cor_matrixOGSIZE[i, j] <- calculate_nonzero_mean_and_wilx2(OGcounts_filt[, i], OGcounts_filt[, j])$mean
    cor_matrixOGSIZETEST[i, j] <- calculate_nonzero_mean_and_wilx2(OGcounts_filt[, i], OGcounts_filt[, j])$pval
    cor_matrixOGSIZECOR[i, j] <- calculate_nonzero_mean_and_wilx2(OGcounts_filt[, i], OGcounts_filt[, j])$cor

  }
}

cor_matrixOGSIZE[lower.tri(cor_matrixOGSIZE)] <- -t(cor_matrixOGSIZE)[lower.tri(cor_matrixOGSIZE)]
cor_matrixOGSIZECOR[lower.tri(cor_matrixOGSIZECOR)] <- t(cor_matrixOGSIZECOR)[lower.tri(cor_matrixOGSIZECOR)]

colnames(cor_matrixOGSIZE) <- sp_real_order
rownames(cor_matrixOGSIZE) <- sp_real_order
colnames(cor_matrixOGSIZECOR) <- sp_real_order
rownames(cor_matrixOGSIZECOR) <- sp_real_order



cor_matrixOGSIZECOR
#C. brenneri C. briggsae C. elegans C. inopinata C. nigoni C. remanei C. tropicalis
#C. brenneri            NA   0.3718833 0.39600599   0.18947921 0.3229793  0.3171207     0.3355575
#C. briggsae     0.3718833          NA 0.53690828   0.20784201 0.5280253  0.4113308     0.5131950
#C. elegans      0.3960060   0.5369083         NA   0.06881321 0.3634646  0.4015517     0.3183854
#C. inopinata    0.1894792   0.2078420 0.06881321           NA 0.1257538  0.2018342     0.1273928
#C. nigoni       0.3229793   0.5280253 0.36346462   0.12575380        NA  0.3267159     0.3458613
#C. remanei      0.3171207   0.4113308 0.40155173   0.20183418 0.3267159         NA     0.3486714
#C. tropicalis   0.3355575   0.5131950 0.31838536   0.12739279 0.3458613  0.3486714            NA

p_values_vectorOG <- as.vector(cor_matrixOGSIZETEST)

# correction
adjusted_p_valuesOG <- p.adjust(p_values_vectorOG, method = "bonferroni")

# back to a matrix
adjOGSIZETEST <- matrix(adjusted_p_valuesOG, nrow = nrow(cor_matrixOGSIZETEST), ncol = ncol(cor_matrixOGSIZETEST))

adjOGSIZETEST[lower.tri(adjOGSIZETEST)] <- t(adjOGSIZETEST)[lower.tri(adjOGSIZETEST)]

colnames(adjOGSIZETEST) <- sp_real_order
rownames(adjOGSIZETEST) <- sp_real_order


colnames(adjOGSIZETEST) <- sp_real_order
rownames(adjOGSIZETEST) <- sp_real_order

adjOGSIZETEST
#C. brenneri  C. briggsae   C. elegans C. inopinata    C. nigoni   C. remanei C. tropicalis
#C. brenneri             NA 6.390706e-05 1.000000e+00 5.004983e-19 4.235080e-02 1.000000e+00  1.918971e-13
#C. briggsae   6.390706e-05           NA 1.205900e-05 2.168523e-09 9.276274e-79 1.178369e-16  1.000000e+00
#C. elegans    1.000000e+00 1.205900e-05           NA 8.660188e-18 3.217510e-02 1.047892e-01  8.184427e-05
#C. inopinata  5.004983e-19 2.168523e-09 8.660188e-18           NA 2.585384e-34 2.874329e-29  9.117562e-08
#C. nigoni     4.235080e-02 9.276274e-79 3.217510e-02 2.585384e-34           NA 1.000000e+00  7.850098e-23
#C. remanei    1.000000e+00 1.178369e-16 1.047892e-01 2.874329e-29 1.000000e+00           NA  5.432028e-17
#C. tropicalis 1.918971e-13 1.000000e+00 8.184427e-05 9.117562e-08 7.850098e-23 5.432028e-17            NA

#I gave up on how to align the title to the left in a normal way, so here we go...
#col.lim = c(-150, 150),
corrplot(cor_matrixOGSIZE[sp_phylo_order, sp_phylo_order], insig = "label_sig",
         sig.level = c(.0001, .001, .01), mar=c(0,0,3,0), title="A. Mean difference in sizes of shared OGs                                                                                  ",  is.corr = FALSE, col.lim = c(-0.3, 0.3), p.mat=adjOGSIZETEST[sp_phylo_order, sp_phylo_order],pch.cex = 0.8, diag = FALSE, font = 3, tl.col = 'black',cl.ratio = 0.25,col=colorRampPalette(c("#633f01","white","#016163"))(6))

#

grab_grob <- function(){
  grid.echo()
  grid.grab()
}


cor_plotOGSIZE <- grab_grob()
cor_plotOGSIZE <- as.ggplot(cor_plotOGSIZE)
cor_plotOGSIZE  <- cor_plotOGSIZE  +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

cor_matrixOGSIZE
#C. brenneri C. briggsae   C. elegans C. inopinata    C. nigoni    C. remanei C. tropicalis
#C. brenneri              NA 0.070785660  0.004132231  0.062969925 -0.026182751  0.0007765784   0.086202912
#C. briggsae   -0.0707856598          NA -0.048482605 -0.011733505 -0.228491903 -0.0901645447  -0.004359591
#C. elegans    -0.0041322314 0.048482605           NA  0.005005247 -0.030762271 -0.0206637445   0.038058252
#C. inopinata  -0.0629699248 0.011733505 -0.005005247           NA -0.008348049 -0.0644997889  -0.005965203
#C. nigoni      0.0261827511 0.228491903  0.030762271  0.008348049           NA  0.0202614379   0.116477483
#C. remanei    -0.0007765784 0.090164545  0.020663745  0.064499789 -0.020261438            NA   0.082315903
#C. tropicalis -0.0862029119 0.004359591 -0.038058252  0.005965203 -0.116477483 -0.0823159026            NA



ratiosOG <- (cor_matrixOGSUM2 / t(cor_matrixOGSUM2) -1 )* 100
ratiosOG
#C. inopinata C. elegans C. tropicalis C. brenneri C. briggsae  C. nigoni  C. remanei
#C. inopinata     0.0000000 -0.4195994    -0.5032502 -5.18575851   1.0167131  -0.676071 -5.24617181
#C. elegans       0.4213674  0.0000000     3.1925984 -0.33350436   4.1324921  -2.452218 -1.63204748
#C. tropicalis    0.5057956 -3.0938250     0.0000000 -6.52973099   0.3645735  -8.821712 -6.25145790
#C. brenneri      5.4693878  0.3346203     6.9858909  0.00000000   5.8894460  -1.982970  0.05808213
#C. briggsae     -1.0064801 -3.9684944    -0.3632492 -5.56188193   0.0000000 -15.792235 -6.95617209
#C. nigoni        0.6806728  2.5138632     9.6752326  2.02308769  18.7538941   0.000000  1.53128430
#C. remanei       5.5366331  1.6591252     6.6683255 -0.05804841   7.4762316  -1.508190  0.00000000



FIGS72 <- cor_plotOGSIZE + cor_plotOGPER

tiff(filename = "FigS7.2.tiff",width =19.6,height = 10.6, res=300, units = "in")  #Figure S5
FIGS72
dev.off()


#### intron sizes
for (i in 1:(ncol(intrones_SIZE) - 1)) {
  for (j in (i + 1):ncol(intrones_SIZE)) {
    cor_matrixSIZE[i, j] <- calculate_nonzero_mean_and_wilx(intrones_SIZE[, i], intrones_SIZE[, j])$mean
    cor_matrixSIZETEST[i, j] <- calculate_nonzero_mean_and_wilx(intrones_SIZE[, i], intrones_SIZE[, j])$pval
    cor_matrixSIZECOR[i, j] <- calculate_nonzero_mean_and_wilx(intrones_SIZE[, i], intrones_SIZE[, j])$cor

  }
}


cor_matrixSIZE[lower.tri(cor_matrixSIZE)] <- -t(cor_matrixSIZE)[lower.tri(cor_matrixSIZE)]
cor_matrixSIZECOR[lower.tri(cor_matrixSIZECOR)] <- t(cor_matrixSIZECOR)[lower.tri(cor_matrixSIZECOR)]

#sp_real_order <- c("C. brenneri","C. briggsae","C. elegans","C. inopinata","C. nigoni","C. remanei","C. tropicalis")
#sp_phylo_order <- c("C. inopinata","C. elegans","C. tropicalis","C. brenneri","C. briggsae","C. nigoni","C. remanei")


colnames(cor_matrixSIZE) <- sp_real_order
rownames(cor_matrixSIZE) <- sp_real_order
colnames(cor_matrixSIZECOR) <- sp_real_order
rownames(cor_matrixSIZECOR) <- sp_real_order

p_values_vector <- as.vector(cor_matrixSIZETEST)

# correction
adjusted_p_values <- p.adjust(p_values_vector, method = "bonferroni")

# back to a matrix
adjSIZETEST <- matrix(adjusted_p_values, nrow = nrow(cor_matrixSIZETEST), ncol = ncol(cor_matrixSIZETEST))

adjSIZETEST[lower.tri(adjSIZETEST)] <- t(adjSIZETEST)[lower.tri(adjSIZETEST)]

colnames(adjSIZETEST) <- sp_real_order
rownames(adjSIZETEST) <- sp_real_order


colnames(adjSIZETEST) <- sp_real_order
rownames(adjSIZETEST) <- sp_real_order


#I gave up on how to align the title to the left in a normal way, so here we go...
corrplot(cor_matrixSIZE[sp_phylo_order, sp_phylo_order], insig = "label_sig",
         sig.level = c(.0001, .001, .01), mar=c(0,0,1,0), title="A. Intron size difference                                                                                               ", col.lim = c(-150, 150), is.corr = FALSE,  p.mat=adjSIZETEST[sp_phylo_order, sp_phylo_order],pch.cex = 0.8, diag = FALSE, font = 3, tl.col = 'black',cl.ratio = 0.25,col=colorRampPalette(c("#633f01","white","#016163"))(10))



grab_grob <- function(){
  grid.echo()
  grid.grab()
}


cor_plotSIZE <- grab_grob()
cor_plotSIZE <- as.ggplot(cor_plotSIZE)



#### now phases and splicing sites (ss)

calculate_nonzero_cramer <- function(col1, col2) {
  nonzero_rows <- which(col1 != "0" & col2 != "0")
  data <- as.data.frame(cbind(as.character(col1[nonzero_rows]),as.character(col2[nonzero_rows])))
  #data <-  apply(data, 2, as.factor)
  if (length(nonzero_rows) > 0) {
    cramer <- cramerV(table(data$V1, data$V2))
    pval <- chisq.test(table(data$V1, data$V2))$p.value
    result <- list(cramer = cramer, pval = pval)
    return(result)
  } else {
    result <- list(cramer = NA, pval = NA)
    return(result)  # Return NA if there are no common non-zero values
  }
}




intrones_PHASE <- as.matrix(intstrALL[, grep("_phase", names(intstrALL))])
intrones_PHASE <-as.data.frame(intrones_PHASE)

intrones_SS <- as.matrix(intstrALL[, grep("_ss", names(intstrALL))])
intrones_SS <- as.data.frame(intrones_SS)


cor_matrixPHASE <- matrix(NA, nrow = ncol(intrones_PHASE), ncol = ncol(intrones_PHASE))
cor_matrixPHASETEST <- matrix(NA, nrow = ncol(intrones_PHASE), ncol = ncol(intrones_PHASE))
for (i in 1:(ncol(intrones_PHASE) - 1)) {
  for (j in (i + 1):ncol(intrones_PHASE)) {
    cor_matrixPHASE[i, j] <- calculate_nonzero_cramer(intrones_PHASE[, i], intrones_PHASE[, j])$cramer
    cor_matrixPHASETEST[i, j] <- calculate_nonzero_cramer(intrones_PHASE[, i], intrones_PHASE[, j])$pval

  }
}


cor_matrixPHASE[lower.tri(cor_matrixPHASE)] <- t(cor_matrixPHASE)[lower.tri(cor_matrixPHASE)]


cor_matrixSS <- matrix(NA, nrow = ncol(intrones_SS), ncol = ncol(intrones_SS))
cor_matrixSSTEST <- matrix(NA, nrow = ncol(intrones_SS), ncol = ncol(intrones_SS))
for (i in 1:(ncol(intrones_SS) - 1)) {
  for (j in (i + 1):ncol(intrones_SS)) {
    cor_matrixSS[i, j] <- calculate_nonzero_cramer(intrones_SS[, i], intrones_SS[, j])$cramer
    cor_matrixSSTEST[i, j] <- calculate_nonzero_cramer(intrones_SS[, i], intrones_SS[, j])$pval
  }
}

cor_matrixSS[lower.tri(cor_matrixSS)] <- t(cor_matrixSS)[lower.tri(cor_matrixSS)]


## need to correct pvals

p_values_vectorPH <- as.vector(cor_matrixPHASETEST)
# correction
adjusted_p_valuesPH <- p.adjust(p_values_vectorPH, method = "bonferroni")
# back to a matrix
adjTESTPH <- matrix(adjusted_p_valuesPH, nrow = nrow(cor_matrixPHASETEST), ncol = ncol(cor_matrixPHASETEST))
adjTESTPH[lower.tri(adjTESTPH)] <- t(adjTESTPH)[lower.tri(adjTESTPH)]

p_values_vectorSS <- as.vector(cor_matrixSSTEST)
# correction
adjusted_p_valuesSS <- p.adjust(p_values_vectorSS, method = "bonferroni")
# back to a matrix
adjTESTSS <- matrix(adjusted_p_valuesSS, nrow = nrow(cor_matrixSSTEST), ncol = ncol(cor_matrixSSTEST))
adjTESTSS[lower.tri(adjTESTSS)] <- t(adjTESTSS)[lower.tri(adjTESTSS)]


colnames(cor_matrixPHASE) <- sp_real_order
colnames(cor_matrixSS)<- sp_real_order
rownames(cor_matrixPHASE)<- sp_real_order
rownames(cor_matrixSS)<- sp_real_order

colnames(adjTESTSS)<- sp_real_order
colnames(adjTESTPH)<- sp_real_order
rownames(adjTESTSS)<- sp_real_order
rownames(adjTESTPH)<- sp_real_order


corrplot(cor_matrixPHASE[sp_phylo_order,sp_phylo_order],insig = "label_sig",sig.level = c(.0001, .001, .01), mar=c(0,0,1,0), title="B. Intron phase                                                                                                ",  is.corr = TRUE,  col.lim=c(0,1),p.mat=adjTESTPH[sp_phylo_order, sp_phylo_order],pch.cex = 0.8, diag = FALSE, font = 3, tl.col = 'black', col=colorRampPalette(c("#633f01","white","#016163"))(10),cl.ratio = 0.25, type = "upper")



cor_plotPH <- grab_grob()
cor_plotPH <- as.ggplot(cor_plotPH)



corrplot(cor_matrixSS[sp_phylo_order,sp_phylo_order], insig = "label_sig",sig.level = c(.0001, .001, .01), mar=c(0,0,1,0), title="C. Intron splicing sites                                                                                                ",  is.corr = TRUE,  col.lim=c(0,1),p.mat=adjTESTSS[sp_phylo_order, sp_phylo_order],pch.cex = 0.8, diag = FALSE, font = 3, tl.col = 'black', col=colorRampPalette(c("#633f01","white","#016163"))(10),cl.ratio = 0.25, type = "upper")


cor_plotSS <- grab_grob()
cor_plotSS <- as.ggplot(cor_plotSS)





FIGS4 <- (cor_plotSIZE + cor_plotPH + cor_plotSS )
tiff(filename = "FigS4.tiff",width = 25.6,height = 9.6, res=300, units = "in") #Figure S9
FIGS4
dev.off()


####################################################
########## Table 2 #################################
####################################################


### stat difference between all genes lengths/intron lengths/exon number between selfers and outcrossers
### all


### 1-to-1

STAT3REP <- STAT3
STAT3REP$Reproduction <- "Selfing"
STAT3REP[STAT3REP$Species %in% c("C. brenneri","C. inopinata","C. nigoni","C. remanei"),]$Reproduction <- "Outcrossing"

table(STAT3REP$Reproduction)

#Outcrossing     Selfing
#1596474     1140331




   calculate_coh_d <- function(df) {
     df %>% coh_d(value ~ Reproduction) -> coh_result

     comparison_direction <- paste0(unique(coh_result$Reproduction_foc),"_",unique(coh_result$Reproduction_ref))
     wt <- wilcox.test(df[df$Reproduction=="Selfing",]$value, df[df$Reproduction=="Outcrossing",]$value, paired = FALSE)$p.value
     return(data.frame(
       Type = df$Type[1],
       Gene= df$Gene[1],
       comp = comparison_direction,
       coh_d_value = coh_result$coh_d,
       wt_pval=wt
     ))
   }

   chdmaGTSTATS <- STAT3REP %>%
     group_split(Type, Gene) %>%
     map_dfr(~ calculate_coh_d(.)) %>% as.data.frame()


   chdmaGTSTATS[chdmaGTSTATS$comp=="Selfing_Outcrossing",] #everything is slightly  longer/more in selfing

   ### Table 2

   #Type   Gene                comp coh_d_value       wt_pval
   #1    gene length 1-to-1 Selfing_Outcrossing 0.069325446  7.278668e-09
   #3    gene length    all Selfing_Outcrossing 0.112834362 1.270538e-156
   #5    exon length 1-to-1 Selfing_Outcrossing 0.056625180 1.272544e-115
   #7    exon length    all Selfing_Outcrossing 0.007340664  5.096278e-81
   #9  intron length 1-to-1 Selfing_Outcrossing 0.009077853 1.317530e-154
   #11 intron length    all Selfing_Outcrossing 0.005825256  0.000000e+00
   #13   exon number 1-to-1 Selfing_Outcrossing 0.049327178  7.271018e-07
   #15   exon number    all Selfing_Outcrossing 0.156966980 3.614828e-299


   library(lmerTest)


   calculate_coh_d2 <- function(df) {
     df %>% coh_d(value ~ Reproduction) -> coh_result

     comparison_direction <- paste0(unique(coh_result$Reproduction_foc), "_", unique(coh_result$Reproduction_ref))

     mixmod<-summary(lmer(value ~ Reproduction + (1|Species),  data=df))

     mixmod_stat <- mixmod$coefficients[2,4]
     mixmod_pval <- mixmod$coefficients[2,5]


     selfing_mean <- mean(df[df$Reproduction == "Selfing", ]$value)
     selfing_sd <- sd(df[df$Reproduction == "Selfing", ]$value)

     outcrossing_mean <- mean(df[df$Reproduction == "Outcrossing", ]$value)
     outcrossing_sd <- sd(df[df$Reproduction == "Outcrossing", ]$value)

     return(data.frame(
       Type = df$Type[1],
       Gene = df$Gene[1],
       comp = comparison_direction,
       coh_d_value = coh_result$coh_d,
       mixmod_stat=mixmod_stat,
       mixmod_pval=mixmod_pval,
       selfing_mean = selfing_mean,
       selfing_sd = selfing_sd,
       outcrossing_mean = outcrossing_mean,
       outcrossing_sd = outcrossing_sd,
       row.names = NULL
     ))
   }


   chdmaGTSTATS2 <- STAT3REP %>%
     group_split(Type, Gene) %>%
     map_dfr(~ calculate_coh_d2(.)) %>% as.data.frame()


   chdmaGTSTATS2[chdmaGTSTATS2$comp=="Selfing_Outcrossing",]
   chdmaGTSTATS2<- chdmaGTSTATS2[chdmaGTSTATS2$comp=="Selfing_Outcrossing",]
   chdmaGTSTATS2$adjpval <- p.adjust(chdmaGTSTATS2$mixmod_pval, method = "bonferroni")
   chdmaGTSTATS2

#  Type   Gene                comp coh_d_value mixmod_stat mixmod_pval selfing_mean  selfing_sd outcrossing_mean outcrossing_sd   adjpval
#   1    gene length 1-to-1 Selfing_Outcrossing 0.069325446   0.5556281  0.60240849  3508.744678 4065.580676      3253.600573    3362.652585 1.0000000
#   3    gene length    all Selfing_Outcrossing 0.112834362   1.0117179  0.35810288  2891.843318 3763.757574      2523.225282    2914.285395 1.0000000
#   5    exon length 1-to-1 Selfing_Outcrossing 0.056625180   0.9606605  0.38084084   218.492615  235.793456       205.474581     225.228476 1.0000000
#   7    exon length    all Selfing_Outcrossing 0.007340664   0.2794799  0.79106126   227.892549  257.277862       225.916684     277.164580 1.0000000
#   9  intron length 1-to-1 Selfing_Outcrossing 0.009077853   0.1551146  0.88279857   306.381496  941.553569       299.049472     685.654392 1.0000000
#   11 intron length    all Selfing_Outcrossing 0.005825256   0.0348734  0.97352998   282.172341  914.221522       277.668941     653.868694 1.0000000
#   13   exon number 1-to-1 Selfing_Outcrossing 0.049327178   1.3159386  0.24529579     7.268650    4.725622         7.041587       4.509204 1.0000000
#   15   exon number    all Selfing_Outcrossing 0.156966980   3.5500181  0.01628827     6.222768    4.378315         5.561903       4.101854 0.1303061


#### test how differenet factor affect the number/size of introns

#intronsize1to1 <- STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="intron length",]
#intronsize1to1$Reproduction <- "Selfing"
#intronsize1to1[intronsize1to1$Species %in% c("C. brenneri","C. inopinata","C. nigoni","C. remanei"),]$Reproduction <- "Outcrossing"


#intronsizeall <- STAT3[STAT3$Gene=="all" & STAT3$Type=="intron length",]
#intronsizeall$Reproduction <- "Selfing"
#intronsizeall[intronsizeall$Species %in% c("C. brenneri","C. inopinata","C. nigoni","C. remanei"),]$Reproduction <- "Outcrossing"


#intronsizeall <- STAT3[STAT3$Gene=="all" & STAT3$Type=="intron length",]

#modelAOVLENG1to1 <- aov(value ~ Reproduction, data = intronsize1to1)


#summary(modelAOVLENG1to1)
#Df    Sum Sq Mean Sq F value Pr(>F)
#Reproduction      1 4.166e+06 4166304   6.387 0.0115 *
#  Residuals    314899 2.054e+11  652351
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#summary(modelAOVLENG1to1)
#Call:
#  aov(formula = value ~ Reproduction, data = intronsize1to1)

#Terms:
#  Reproduction    Residuals
#Sum of Squares       4166304 205424760557
#Deg. of Freedom            1       314899

#Residual standard error: 807.6827
#Estimated effects may be unbalanced

#summary(intronsize1to1[intronsize1to1$Reproduction=="Selfing",]$value)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#1.0     46.0     53.0    306.4    269.0 137839.0
#summary(intronsize1to1[intronsize1to1$Reproduction=="Outcrossing",]$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#29.0    47.0    55.0   299.1   260.0 47775.0
#coh_d(value~ Reproduction, data=intronsize1to1)
#1 Outcrossing      Selfing           0.00908 0.00359
#2 Selfing          Outcrossing      -0.00908 0.00359



#intronsize1to1 %>% coh_d(value~ Reproduction) -> chdREP
## A tibble: 2 × 4
#Reproduction_ref Reproduction_foc    coh_d  coh_se
#<chr>            <chr>               <dbl>   <dbl>
#  1 Outcrossing      Selfing           0.00908 0.00359
#2 Selfing          Outcrossing      -0.00908 0.00359


#modelAOVLENGall <- aov(value ~ Reproduction, data = intronsizeall)

#summary(aov(value ~ Reproduction, data = intronsizeall))
#Df    Sum Sq Mean Sq F value Pr(>F)
#Reproduction      1 3.671e+06 3670524   6.142 0.0132 *
#  Residuals    744624 4.450e+11  597656
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#intronsizeall %>% coh_d(value~ Reproduction)
# A tibble: 2 × 4
#Reproduction_ref Reproduction_foc    coh_d  coh_se
#<chr>            <chr>               <dbl>   <dbl>
#  1 Outcrossing      Selfing           0.00583 0.00235
#2 Selfing          Outcrossing      -0.00583 0.00235



### ok, let's addd random effects:


#library(lme4)

#STAT5 <-STAT3
#STAT5$Reproduction <- "Selfing"
#STAT5[STAT5$Species %in% c("C. brenneri","C. inopinata","C. nigoni","C. remanei"),]$Reproduction <- "Outcrossing"

#summary(glm(value ~ Reproduction + (1|Species), family = gaussian(link = "identity"), data=STAT5[STAT5$Gene=="1-to-1" & STAT5$Type=="intron length",]))






num_cols <- ncol(intrones_SIZE)


  resultIntDif <- data.frame()

  for (i in 1:(num_cols - 1)) {
    for (j in (i + 1):num_cols) {
      col1 <- intrones_SIZE[,i]
      col2 <- intrones_SIZE[,j]
      nonzero_rows <- which(col1 > 0 & col2 > 0)
      data <- as.data.frame(cbind(col1[nonzero_rows],col2[nonzero_rows]))

      if (length(nonzero_rows) > 0) {

        resultIntDif <- rbind(resultIntDif, data.frame(Column1 = colnames(intrones_SIZE)[i], Column2 = colnames(intrones_SIZE)[j], Difference = data$V1 - data$V2))

      }
    }
  }

head(resultIntDif)

resultIntDif$Column1 <- gsub("_size","",resultIntDif$Column1)
resultIntDif$Column2 <- gsub("_size","",resultIntDif$Column2)
resultIntDif$Reproduction <- "Selfing"
resultIntDif$Reproduction2 <- "Selfing"
resultIntDif[resultIntDif$Column1 %in% c("Cbren","Cinop","Crema"),]$Reproduction <- "Outcrossing"
resultIntDif[resultIntDif$Column2 %in% c("Cbren","Cinop","Crema"),]$Reproduction2 <- "Outcrossing"

#TREESP<-c("(Cinop:0.091484,(Celeg:0.134144,((Crema:0.127918,(Cnigo:0.0504089,Cbrig:0.0334806)0.90561:0.117991)0.440453:0.026276,(Ctrop:0.129409,Cbren:0.125876)0.533563:0.0320385)0.636122:0.036656)1:0.091484);")
melted1x <- melted1
colnames(melted1x) <-c("Column1","Column2","phylo_distance")

resultIntDifP <- merge(melted1x,resultIntDif,by=c("Column1","Column2"))

#lets flip some to have only 3 classes
flip_condition <- resultIntDifP$Reproduction == "Outcrossing" & resultIntDifP$Reproduction2 == "Selfing"
resultIntDifP2 <- resultIntDifP

resultIntDifP2$Column1[flip_condition] <- resultIntDifP$Column2[flip_condition]
resultIntDifP2$Column2[flip_condition] <- resultIntDifP$Column1[flip_condition]
resultIntDifP2$Difference[flip_condition] <- -resultIntDifP$Difference[flip_condition]
resultIntDifP2$Reproduction[flip_condition] <- resultIntDifP$Reproduction2[flip_condition]
resultIntDifP2$Reproduction2[flip_condition] <- resultIntDifP$Reproduction[flip_condition]
resultIntDifP2$Pair <- paste0(resultIntDifP2$Reproduction,"_",resultIntDifP2$Reproduction2)


resultIntDifP2 %>% coh_d(Difference~ Pair)
## A tibble: 6 × 4
#Pair_ref                Pair_foc                  coh_d  coh_se
#<chr>                   <chr>                     <dbl>   <dbl>
#  1 Outcrossing_Outcrossing Selfing_Outcrossing     -0.0439 0.00559
#2 Outcrossing_Outcrossing Selfing_Selfing          0.0697 0.00600
#3 Selfing_Outcrossing     Outcrossing_Outcrossing  0.0439 0.00559
#4 Selfing_Outcrossing     Selfing_Selfing          0.117  0.00405
#5 Selfing_Selfing         Outcrossing_Outcrossing -0.0697 0.00600
#6 Selfing_Selfing         Selfing_Outcrossing     -0.117  0.00405

mean_differencesREP <- resultIntDifP2 %>%
  group_by(Reproduction, Reproduction2) %>%
  summarise(mean_difference = mean(Difference, na.rm = TRUE))
#summarise()` has grouped output by 'Reproduction'. You can override using the `.groups` argument.
# A tibble: 3 × 3
# Groups:   Reproduction [2]
#Reproduction Reproduction2 mean_difference
#<chr>        <chr>                   <dbl>
#  1 Outcrossing  Outcrossing              14.3
#2 Selfing      Outcrossing             -18.3
#3 Selfing      Selfing                  71.5


# Fit a linear model with Species as random effects

library(MASS)
library(lme4)

robust_modelREPPHY <- rlm(Difference ~ phylo_distance + Reproduction, data = resultIntDifP2)
robust_residualsREPPHY <- residuals(robust_modelREPPHY)


mixed_modelREPPHY <- lmer(Difference ~ phylo_distance + Reproduction + (1|Column1) + (1|Column2), data = resultIntDifP2)
mixed_modelREPPHYRES <- lmer(Difference ~ phylo_distance + Reproduction + (1|Column1) + (1|Column2), data = resultIntDifP2, weights = 1/abs(robust_residualsREPPHY))



summary(mixed_modelREPPHY)
#Formula: Difference ~ phylo_distance + Reproduction + (1 | Column1) +      (1 | Column2)
#Data: resultIntDifP2
#
#REML criterion at convergence: 4933767
#
#Scaled residuals:
#Min      1Q  Median      3Q     Max
#-62.396  -0.125  -0.039   0.107 112.100

#Random effects:
#Groups   Name        Variance Std.Dev.
#Column1  (Intercept)   4678    68.39
#Column2  (Intercept)   5721    75.64
#Residual             586639   765.92
#Number of obs: 306062, groups:  Column1, 6; Column2, 6

#Fixed effects:
#Estimate Std. Error         df t value Pr(>|t|)
#(Intercept)             13.704     58.470      7.484   0.234   0.8210
#phylo_distance          66.071     28.803 176067.742   2.294   0.0218 *
#ReproductionSelfing    -39.198     59.431      4.018  -0.660   0.5454
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#(Intr) phyl_d
#phylo_dstnc -0.174
#RprdctnSlfn -0.682  0.018
summary(mixed_modelREPPHYRES)

#Formula: Difference ~ phylo_distance + Reproduction + (1 | Column1) +      (1 | Column2)
#Data: resultIntDifP2
#Weights: 1/abs(robust_residualsREPPHY)

#REML criterion at convergence: 3474966

#Scaled residuals:
#  Min       1Q   Median       3Q      Max
#-14.6800  -0.2420  -0.0382   0.2786  19.6939

#Random effects:
#  Groups   Name        Variance Std.Dev.
#Column1  (Intercept)   0.5023  0.7087
#Column2  (Intercept)   0.4314  0.6568
#Residual             221.4737 14.8820
#Number of obs: 306062, groups:  Column1, 6; Column2, 6

#Fixed effects:
#  Estimate Std. Error       df t value Pr(>|t|)
#(Intercept)           1.1234     0.6821  11.4991   1.647    0.127
#phylo_distance       -1.2420     0.9707 684.6662  -1.279    0.201
#ReproductionSelfing  -0.4827     0.6401   4.3146  -0.754    0.490
#
#Correlation of Fixed Effects:
#  (Intr) phyl_d
#phylo_dstnc -0.487
#RprdctnSlfn -0.669  0.053

AIC(mixed_modelREPPHY) < AIC(mixed_modelREPPHYRES) #FALSE
BIC(mixed_modelREPPHY) < BIC(mixed_modelREPPHYRES) #FALSE




###ok . let's try to treat species as a random effeect
#https://people.bath.ac.uk/jjf23/mixchange/rbd.html



  ############ ok. now number of introns


intronnum1to1 <- STAT3[STAT3$Gene=="1-to-1" & STAT3$Type=="exon number",]
intronnum1to1$Reproduction <- "Selfing"
intronnum1to1[intronnum1to1$Species %in% c("C. brenneri","C. inopinata","C. nigoni","C. remanei"),]$Reproduction <- "Outcrossing"


intronnumall <- STAT3[STAT3$Gene=="all" & STAT3$Type=="exon number",]
intronnumall$Reproduction <- "Selfing"
intronnumall[intronnumall$Species %in% c("C. brenneri","C. inopinata","C. nigoni","C. remanei"),]$Reproduction <- "Outcrossing"





coh_d(value~ Reproduction, data=intronnum1to1)
## A tibble: 2 × 4
#Reproduction_ref Reproduction_foc   coh_d  coh_se
#<chr>            <chr>              <dbl>   <dbl>
#  1 Outcrossing      Selfing           0.0493 0.00892
#2 Selfing          Outcrossing      -0.0493 0.00892


coh_d(value~ Reproduction, data=intronnumall)
# A tibble: 2 × 4
#Reproduction_ref Reproduction_foc  coh_d  coh_se
#<chr>            <chr>             <dbl>   <dbl>
#  1 Outcrossing      Selfing           0.157 0.00524
#2 Selfing          Outcrossing      -0.157 0.00524


modelAOVNUM1to1 <- aov(value ~ Reproduction, data = intronnum1to1)
summary(modelAOVNUM1to1)
#Df  Sum Sq Mean Sq F value   Pr(>F)
#Reproduction     1     648   647.7   30.57 3.24e-08 ***
#  Residuals    51294 1086892    21.2
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



modelAOVNUMall <- aov(value ~ Reproduction, data = intronnumall)
summary(modelAOVNUMall)
#Df  Sum Sq Mean Sq F value Pr(>F)
#Reproduction      1   15977   15977   901.3 <2e-16 ***
#  Residuals    154619 2740771      18
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



summary(intronnum1to1[intronnum1to1$Reproduction=="Selfing",]$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   4.000   6.000   7.269   9.000  64.000

summary(intronnum1to1[intronnum1to1$Reproduction=="Outcrossing",]$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   4.000   6.000   7.042   9.000  61.000

summary(intronnumall[intronnumall$Reproduction=="Selfing",]$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   3.000   5.000   6.223   8.000  67.000

summary(intronnumall[intronnumall$Reproduction=="Outcrossing",]$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   3.000   5.000   5.562   7.000  65.000

sd(intronnumall[intronnumall$Reproduction=="Selfing",]$value)
#[1] 4.378315

sd(intronnumall[intronnumall$Reproduction=="Outcrossing",]$value)
#[1] 4.101854

sd(intronnum1to1[intronnum1to1$Reproduction=="Selfing",]$value)
#[1] 4.725622

sd(intronnum1to1[intronnum1to1$Reproduction=="Outcrossing",]$value)
#[1] 4.509204



##############################################
################## Fig8 ######################
##############################################


intrPHASE <- as.matrix(intstrALL[, grep("_phase", names(intstrALL))])


#intrones_PHASE <-as.data.frame(intrones_PHASE)
#and intrones_SS


#remove singletons
intrPHASECOM<-intrones_PHASE[which(rowSums(intrones_PHASE != "0") > 1),]
intrPHASEUNI<-intrones_PHASE[which(rowSums(intrones_PHASE != "0") == 1),]

intrSSCOM<-intrones_SS[which(rowSums(intrones_SS != "0") > 1),]
intrSSUNI<-intrones_SS[which(rowSums(intrones_SS != "0") == 1),]
intrSSUNI<-data.frame(subset(intrSSUNI,!apply(intrSSUNI, 1, function(row) any(grepl("N", row)))))
intrSSCOM<-data.frame(subset(intrSSCOM,!apply(intrSSCOM, 1, function(row) any(grepl("N", row)))))

#subset(df, !apply(df, 1, function(row) any(grepl("N", row))))

intrPHASECOM_sum <- intrPHASECOM %>%
  mutate_all(as.character) %>%
  gather(column, value) %>%
  filter(value != "0") %>%
  group_by(column, value) %>%
  summarise(count = n()) %>% as.data.frame()

intrPHASEUNI_sum <- intrPHASEUNI %>%
  mutate_all(as.character) %>%
  gather(column, value) %>%
  filter(value != "0") %>%
  group_by(column, value) %>%
  summarise(count = n()) %>% as.data.frame()


intrSSCOM_sum <- intrSSCOM %>%
  mutate_all(as.character) %>%
  gather(column, value) %>%
  filter(value != "0") %>%
  group_by(column, value) %>%
  summarise(count = n()) %>% as.data.frame()


intrSSUNI_sum <- intrSSUNI %>%
  mutate_all(as.character) %>%
  gather(column, value) %>%
  filter(value != "0") %>%
  group_by(column, value) %>%
  summarise(count = n()) %>% as.data.frame()



intrPHASECOM_sum$Type <- "Common"
intrPHASEUNI_sum$Type <- "Unique"

intrSSCOM_sum$Type <- "Common"
intrSSUNI_sum$Type <- "Unique"


intrPHASECOM_sum$Class <- gsub("p","Phase ", intrPHASECOM_sum$value)
intrPHASEUNI_sum$Class <- gsub("p","Phase ", intrPHASEUNI_sum$value)

intrSSCOM_sum$Class <- "Other"
intrSSCOM_sum[intrSSCOM_sum$value=="GT-AG",]$Class <- "GT-AG"
intrSSCOM_sum[intrSSCOM_sum$value=="GC-AG",]$Class <- "GC-AG"


intrSSUNI_sum$Class <- "Other"
intrSSUNI_sum[intrSSUNI_sum$value=="GT-AG",]$Class <- "GT-AG"
intrSSUNI_sum[intrSSUNI_sum$value=="GC-AG",]$Class <- "GC-AG"


intrPHASECOMBO <- rbind(intrPHASECOM_sum,intrPHASEUNI_sum)
intrSSCOMBO <- rbind(intrSSCOM_sum,intrSSUNI_sum)


#       column value count   Type   Class
#1 Cbren_phase    p0 14989 Common Phase 0
#2 Cbren_phase    p1  6978 Common Phase 1
#3 Cbren_phase    p2  6954 Common Phase 2




#group_by(Type, column) %>%
#  mutate(percentage = count / sum(count) * 100)


intrPHASECOMBOP <-  intrPHASECOMBO  %>%
  group_by(Type,column) %>%
  mutate(percentage = (count / sum(count)) * 100) %>% as.data.frame()

intrSSCOMBOP <-  intrSSCOMBO  %>%
  group_by(Type,column) %>%
  mutate(percentage = (count / sum(count)) * 100) %>% as.data.frame()


intrSSCOMBOP$column <- gsub("_ss","",intrSSCOMBOP$column)
intrPHASECOMBOP$column <- gsub("_phase","",intrPHASECOMBOP$column)


sublist2<-list("Cbren"="C. brenneri", "Cbrig"="C. briggsae", "Cnigo"="C. nigoni", "Ctrop"="C. tropicalis",  "Cinop"="C. inopinata", "Crema"="C. remanei", "Celeg"="C. elegans")

intrSSCOMBOP$column <- Reduce(function(col, pattern) gsub(pattern, sublist2[[pattern]], col), names(sublist2), init = intrSSCOMBOP$column)
intrPHASECOMBOP$column <- Reduce(function(col, pattern) gsub(pattern, sublist2[[pattern]], col), names(sublist2), init = intrPHASECOMBOP$column)


intrPHASECOMBOP$column <- factor(intrPHASECOMBOP$column, levels = rev(c("C. brenneri","C. remanei","C. nigoni","C. inopinata","C. briggsae","C. tropicalis","C. elegans")))
intrSSCOMBOP$column <- factor(intrSSCOMBOP$column, levels = rev(c("C. brenneri","C. remanei","C. nigoni","C. inopinata","C. briggsae","C. tropicalis","C. elegans")))



intrPHASECOMBOP$Class2 <- "Phase 0"
intrPHASECOMBOP[intrPHASECOMBOP$Class !="Phase 0",]$Class2 <- "not Phase 0"



intrPHASECOMBOP2 <- intrPHASECOMBOP %>%
  group_by(column, Type, Class2) %>%
  summarise(sum_count = sum(count))


perform_fisher_test <- function(count_table) {
  fisher_result <- fisher.test(matrix(count_table, ncol = 2))
  return(fisher_result$p.value)
}






#calculate_odds_ratio <- function(a, b, c, d) {
#  odds_ratio <- (a / b) / (c / d)
#  return(odds_ratio)
#}





intrPHASECOMBOP2FISHER <- intrPHASECOMBOP2  %>%
  group_by(column) %>%
  summarise(
    count_not_Phase0_Common = sum(sum_count[Class2 == "not Phase 0" & Type == "Common"]),
    count_Phase0_Common = sum(sum_count[Class2 == "Phase 0" & Type == "Common"]),
    count_not_Phase0_Unique = sum(sum_count[Class2 == "not Phase 0" & Type == "Unique"]),
    count_Phase0_Unique = sum(sum_count[Class2 == "Phase 0" & Type == "Unique"]),
    p_value = perform_fisher_test(c(
      count_Phase0_Common,
      count_not_Phase0_Common,
      count_Phase0_Unique,
      count_not_Phase0_Unique
    ))
  )

intrPHASECOMBOP2FISHER <-as.data.frame(intrPHASECOMBOP2FISHER)
intrPHASECOMBOP2FISHER$padj <- p.adjust(intrPHASECOMBOP2FISHER$p_value, method = "bonferroni")





intrPHASECOMBOP2<- merge(intrPHASECOMBOP, intrPHASECOMBOP2FISHER[,c(1,7)],by="column")
intrPHASECOMBOP2[intrPHASECOMBOP2$Type=="Unique",]$padj <- NA

phaseper2 <- ggplot(intrPHASECOMBOP2, aes(x = percentage, y=Type, ordered = FALSE, fill=Class)) + theme_bw(base_size = 10) + labs(
  title = "A",
  x = "%",
  y = "") + theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    plot.title = element_text(hjust = 0,face="bold")
  ) + theme(legend.position = "top") + theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = rel(1))
  )  + geom_bar(stat = "identity", position = "stack") +  scale_fill_manual(values = c("#3f918a","#ff7f00","grey20"), name="") +  theme(strip.text = element_text(size = 10))  + theme(strip.text = element_text(size = 10), panel.border = element_blank(), axis.line = element_line(colour = "grey40")) +
  facet_grid( column ~ .) +
  theme(strip.text.y = element_text(face = "italic",angle = 0,hjust=0)) +
  geom_text(data=intrPHASECOMBOP2[intrPHASECOMBOP2$Type == "Common",],aes(x = 54, label = ifelse(padj < 0.05, "]","")),  vjust = -0.1,  # Adjust the position between bars
                                                                                    hjust = -0.40,
                                                                                    size = 6,
                                                                                    color = "grey25") +
  geom_text(
    aes(x = 57, label = ifelse(padj < 0.0001, "***",
                               ifelse(padj < 0.001, "**",
                                      ifelse(padj < 0.01, "*",
                                             ifelse(padj < 0.05, ".", ""))))),
    vjust = 0.1,  # Adjust the position between bars
    hjust = 0,
    size = 4,
    color = "grey25"
  )



#### now ss
intrSSCOMBOP$Class2 <- "GT-AG"
intrSSCOMBOP[intrSSCOMBOP$Class !="GT-AG",]$Class2 <- "not GT-AG"



intrSSCOMBOP2 <- intrSSCOMBOP %>%
  group_by(column, Type, Class2) %>%
  summarise(sum_count = sum(count))





intrSSCOMBOP2FISHER <- intrSSCOMBOP2 %>%
  group_by(column) %>%
  summarise(
    count_not_GT_AG_Common = sum(sum_count[Class2 == "not GT-AG" & Type == "Common"]),
    count_GT_AG_Common = sum(sum_count[Class2 == "GT-AG" & Type == "Common"]),
    count_not_GT_AG_Unique = sum(sum_count[Class2 == "not GT-AG" & Type == "Unique"]),
    count_GT_AG_Unique = sum(sum_count[Class2 == "GT-AG" & Type == "Unique"]),
    p_value = perform_fisher_test(c(
      count_GT_AG_Common,
      count_not_GT_AG_Common,
      count_GT_AG_Unique,
      count_not_GT_AG_Unique
    ))
  )


intrSSCOMBOP2FISHER <-as.data.frame(intrSSCOMBOP2FISHER)


intrSSCOMBOP2FISHER$padj <- p.adjust(intrSSCOMBOP2FISHER$p_value, method = "bonferroni")

intrSSCOMBOP2<- merge(intrSSCOMBOP, intrSSCOMBOP2FISHER[,c(1,7)],by="column")
intrSSCOMBOP2[intrSSCOMBOP2$Type=="Unique",]$padj <- NA




#ssper <-

ssper2 <-   ggplot(intrSSCOMBOP2[intrSSCOMBOP2$Class !="GT-AG",], aes(x = percentage, y=Type, ordered = FALSE, fill=Class)) + theme_bw(base_size = 10) + labs(
  title = "B",
  x = "%",
  y = "") + theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 0.5
    ),
    plot.title = element_text(hjust = 0,face="bold")
  ) + theme(legend.position = "top") + theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = rel(1))
  )  + geom_bar(stat = "identity", position = "stack") +  scale_fill_manual(values = c("#3f918a","#ff7f00","grey95"), name="") +  theme(strip.text = element_text(size = 10))  + theme(strip.text = element_text(size = 10), panel.border = element_blank(), axis.line = element_line(colour = "grey40")) +
  facet_grid( column ~ .) +
  theme(strip.text.y = element_text(face = "italic",angle = 0,hjust=0)) +
  geom_text(data=intrSSCOMBOP2[intrSSCOMBOP2$Type == "Common" & intrSSCOMBOP2$Class !="GT-AG",],aes(x = 5.6, label = ifelse(padj < 0.05, "]","")),  vjust = -0.1,  # Adjust the position between bars
            hjust = -0.40,
            size = 6,
            color = "grey25") +
  geom_text(
    aes(x = 5.6, label = ifelse(padj < 0.0001, "***",
                               ifelse(padj < 0.001, "**",
                                      ifelse(padj < 0.01, "*",
                                             ifelse(padj < 0.05, ".", ""))))),
    vjust = 0.1,  # Adjust the position between bars
    hjust = -0.5,
    size = 4,
    color = "grey25"
  ) + xlim(0,5.8)


FIG6.2 <- phaseper2 + ssper2
tiff(filename = "Fig6.2.tiff",width = 10.4,height = 4.2, res=300, units = "in")  ###Figure 8
plot(FIG6.2)
dev.off()



##############################################
########### Fig5,  FigS3, FigS4 ##############
##############################################


SYN<-read.csv("GENESPACE_exon_intron_frac_count_Cbren_vs_othersFIX.txt",header=TRUE, sep="\t")

SYN2<-SYN[SYN$genome1=="Cbren" & SYN$genome2!="Cbren",]
table(SYN2$genome2)




SYN2$Species <- SYN2$genome2

SYN2$Species <-gsub("Cbren","C. brenneri", SYN2$Species)
SYN2$Species <-gsub("Celeg","C. elegans", SYN2$Species)
SYN2$Species <-gsub("Crema","C. remanei", SYN2$Species)
SYN2$Species <-gsub("Cbrig","C. briggsae",SYN2$Species)
SYN2$Species <-gsub("Cinop","C. inopinata",SYN2$Species)
SYN2$Species <-gsub("Ctrop", "C. tropicalis",SYN2$Species)
SYN2$Species <-gsub("Cnigo", "C. nigoni",SYN2$Species)
#"Cbren" "Cbrig" "Celeg" "Cinop" "Cnigo" "Crema" "Ctrop"

SYN2$Species <- factor(SYN2$Species, levels = c("C. brenneri","C. remanei","C. nigoni","C. briggsae","C. tropicalis","C. elegans","C. inopinata"))

#pseudo count for fully covered areas
SYN2[SYN2$uncovered_region_genome1<0,]$uncovered_region_genome1 <- 1
SYN2[SYN2$uncovered_region_genome2<0,]$uncovered_region_genome2 <- 1
SYN3<- SYN2[!(SYN2$chr1=="V" & SYN2$chr2=="CM008509.1"), ]
SYN3 <- SYN3[!(SYN3$chr1=="IV" & SYN3$chr2=="II"),]

melted1[melted1$Var1=="Cbren",]
#Var1  Var2     value
#7  Cbren Cinop 0.3775385
#14 Cbren Celeg 0.3287145
#21 Cbren Crema 0.3121085
#28 Cbren Cnigo 0.3525904
#35 Cbren Cbrig 0.3356621
#42 Cbren Ctrop 0.2552850

table(SYN3$genome2)

#-   +
#  Cbren   0   0
#Cbrig  65  65
#Celeg 105 116
#Cinop  87  92
#Cnigo  62  70
#Crema  70  61
#Ctrop  40  62




SYN3 %>%
  group_by(genome2) %>%
  filter(endBp1 - startBp1 == max(endBp1 - startBp1, na.rm = TRUE)) %>%
  summarise(max_diff = max(endBp1 - startBp1, na.rm = TRUE),
            corresponding_diff = max(endBp2 - startBp2, na.rm = TRUE))
#genome2 max_diff corresponding_diff
#<fct>      <int>              <int>
#  1 Cbrig   11574891            9138501
#2 Celeg    5931714            3871936
#3 Cinop    4212007             961204
#4 Cnigo   14636720           12767482
#5 Crema   14608066           17085474
#6 Ctrop   23304476           12757971



SYN3 %>%
  group_by(genome2) %>%
  summarise(mean_diff1 = mean(endBp1 - startBp1, na.rm = TRUE)/100000,
            sd_diff1 = sd(endBp1 - startBp1, na.rm = TRUE)/100000,
            mean_diff2 = mean(endBp2 - startBp2, na.rm = TRUE)/100000,
            sd_diff2 = sd(endBp2 - startBp2, na.rm = TRUE)/100000)
## A tibble: 6 × 5
#genome2 mean_diff1 sd_diff1 mean_diff2 sd_diff2
#<fct>        <dbl>    <dbl>      <dbl>    <dbl>
#  1 Cbrig         8.65    16.0        7.89    11.8
#2 Celeg         5.00     7.96       4.23     5.64
#3 Cinop         5.85     7.27       6.37     6.56
#4 Cnigo         8.45    17.0        8.59    15.4
#5 Crema         8.34    15.6        8.72    16.9
#6 Ctrop        11.2     25.9        7.51    15.4


par(mar = c(5.1, 4.1, 4.1, 2.1) + 0.1)




SYN3 %>%
  group_by(genome2,orient) %>%
  summarise(mean_diff1 = mean(endBp1 - startBp1, na.rm = TRUE)/100000,
            sd_diff1 = sd(endBp1 - startBp1, na.rm = TRUE)/100000,
            mean_diff2 = mean(endBp2 - startBp2, na.rm = TRUE)/100000,
            sd_diff2 = sd(endBp2 - startBp2, na.rm = TRUE)/100000)

### Groups:   genome2 [6]
#genome2 orient mean_diff1 sd_diff1 mean_diff2 sd_diff2
#<fct>   <fct>       <dbl>    <dbl>      <dbl>    <dbl>
#  1 Cbrig   -            5.09     4.48       5.34     4.82
#2 Cbrig   +           12.2     21.7       10.4     15.6
#3 Celeg   -            4.75     8.42       3.83     5.37
#4 Celeg   +            5.23     7.54       4.58     5.87
#5 Cinop   -            7.13     9.11       6.95     6.94
#6 Cinop   +            4.63     4.68       5.83     6.18
#7 Cnigo   -            5.50     5.65       6.05     5.88
#8 Cnigo   +           11.1     22.4       10.8     20.2
#9 Crema   -            6.37     7.68       6.76     8.34
#10 Crema   +           10.6     21.2       11.0     23.0
#11 Ctrop   -            6.76     9.56       4.33     5.07
#12 Ctrop   +           14.1     32.2        9.56    19.1





RATIOS <-c()


ratios_to_calculate <- c("exon_length", "intron_length", "uncovered_region", "genecount","exoncount", "introncount")



for (ratio in ratios_to_calculate){

  TMP <-SYN2[,c(1:9,22)]
  TMP$rat_type <- paste0("ratio_", ratio)
  TMP$value <- SYN2[,paste0(ratio, "_genome2")]  / SYN2[,paste0(ratio, "_genome1")]

  RATIOS <-rbind(RATIOS,TMP)

}
TMP <-SYN2[,c(1:9,22)]
TMP$rat_type <- "ratio_size"
TMP$value <- (SYN2$endBp2 - SYN2$startBp2)  / (SYN2$endBp1 - SYN2$startBp1)
RATIOS <-rbind(RATIOS,TMP)

RATIOS$Type <-RATIOS$rat_type
RATIOS$Type <- gsub("ratio_","",RATIOS$Type)
RATIOS$Type <-gsub("count"," count", RATIOS$Type)
RATIOS$Type <-gsub("_length"," length", RATIOS$Type)
RATIOS$Type <-gsub("uncovered_region", "intergenic",RATIOS$Type)
RATIOS$Type <-gsub("size","block size", RATIOS$Type)


RATIOS$Type <- factor(RATIOS$Type, levels=c("block size","gene count","exon count","intron count","exon length","intron length","intergenic"))

#remove 2 strange blocks that are not at the same chromosomes
RATIOS <- RATIOS[!(RATIOS$chr1=="V" & RATIOS$chr2=="CM008509.1"), ]
RATIOS <- RATIOS[!(RATIOS$chr1=="IV" & RATIOS$chr2=="II"),]


ratioexinCOUNT<-ggplot(RATIOS[RATIOS$rat_type %in% c("ratio_exoncount","ratio_introncount"),], aes(x = startBp1, y = log2(value), ordered = FALSE, col=Type,fill=Type)) + facet_grid(Species ~ chr1, scales="free_x") + theme_bw(base_size = 12) + labs(title = "A", x = expression(~ italic("C. brenneri") ~ " genome position (Mb)"), y = expression("Ratio " ~ italic("C.sp") ~"/" ~italic("C. brenneri"))) + theme(
  axis.text.x = element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.title = element_text(face = "bold", hjust = 0)
) + theme(legend.position = "top") + scale_colour_manual(values = c("#FF7C01", "#3bb3ad", "#696866"), name="") + scale_fill_manual(values = c("#FF7C01", "#3bb3ad", "#696866"),name="")  + theme(
  strip.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))
)  +  theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12)) + scale_x_continuous(limits = c(0,NA),
                                                                                                                      labels = function(x)
                                                                                                                        x / 1000000
) + coord_cartesian(clip = "off")  + geom_hline(yintercept = log2(1), linetype = "dotted", linewidth=0.9, color = "black") + geom_rect(aes(xmin = startBp1, xmax = endBp1, ymin=log2(value)-0.2, ymax=log2(value)+0.2), linewidth=1.35, alpha = 0.75, color='NA') + theme(strip.text.y = element_text(face = "italic"))
###
ratioexinCOUNT

##


ratioexinLENGTH<-ggplot(RATIOS[RATIOS$rat_type %in% c("ratio_exon_length","ratio_intron_length","ratio_uncovered_region") & RATIOS$value>0.00001 & RATIOS$value<14,], aes(x = startBp1, y = value, ordered = FALSE, col=Type,fill=Type)) + facet_grid(Species ~ chr1, scales="free_x") + theme_bw(base_size = 12) + labs(title = "B", x = expression(~ italic("C. brenneri") ~ " genome position (Mb)"), y = expression("Ratio " ~ italic("C.sp") ~"/" ~italic("C. brenneri"))) + theme(
  axis.text.x = element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.title = element_text(face = "bold", hjust = 0)
) + theme(legend.position = "top") + scale_colour_manual(values = c("#FF7C01", "#3bb3ad", "#696866"), name="") + scale_fill_manual(values = c("#FF7C01", "#3bb3ad", "#696866"),name="")  + theme(
  strip.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))
)  +  theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12)) + scale_x_continuous(limits = c(0,NA),
                                                                                                                      labels = function(x)
                                                                                                                        x / 1000000
) + coord_cartesian(clip = "off")  + geom_hline(yintercept = log2(1), linetype = "dotted", linewidth=0.9, color = "black") + geom_rect(aes(xmin = startBp1, xmax = endBp1, ymin=log2(value)-0.25, ymax=log2(value)+0.25), linewidth=1.35, alpha = 0.75, color='NA') + theme(strip.text.y = element_text(face = "italic")) + ylim(-5,5)




ratioexinGENE<-    ggplot(RATIOS[RATIOS$rat_type %in% c("ratio_size","ratio_genecount"),], aes(x = startBp1, y = value, ordered = FALSE, col=Type,fill=Type)) + facet_grid(Species ~ chr1, scales="free_x") + theme_bw(base_size = 12) + labs(title = "C", x = expression(~ italic("C. brenneri") ~ " genome position (Mb)"), y = expression("Ratio " ~ italic("C.sp") ~"/" ~italic("C. brenneri"))) + theme(
  axis.text.x = element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.title = element_text(face = "bold", hjust = 0)
) + theme(legend.position = "top") + scale_colour_manual(values = c("#808687","#de7702"), name="") + scale_fill_manual(values = c("#808687", "#de7702"),name="")  + theme(
  strip.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))
)  +  theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12)) + scale_x_continuous(limits = c(0,NA),
                                                                                                                      labels = function(x)
                                                                                                                        x / 1000000
) + coord_cartesian(clip = "off")  + geom_hline(yintercept = log2(1), linetype = "dotted", linewidth=0.9, color = "black") + geom_rect(aes(xmin = startBp1, xmax = endBp1, ymin=log2(value)-0.2, ymax=log2(value)+0.2), linewidth=1.35, alpha = 0.75, color='NA') + theme(strip.text.y = element_text(face = "italic")) + ylim(-4,4)





FIGS5<- (ratioexinCOUNT + ratioexinLENGTH + ratioexinGENE) #Figure S3
tiff(filename = "FigS5.tiff",width = 22.6,height = 10.6, res=300, units = "in")
plot(FIGS5)
dev.off()

#/ratioexinGENE goes to the main figure




###### can we classify selfers and outcrossers and if yes, what the main difference?
library(glmnet)
library(xgboost)
library(tidyr)
library(dplyr)
library(caret)

melted1 -> phylo
colnames(phylo) <-c("genome1","genome2","phylo_dist")



RATIOS2 <- RATIOS[, -which(names(RATIOS) == "Type")]

RATIOS_WIDE <- pivot_wider(
  data = RATIOS2,
  id_cols = c("genome1", "genome2", "chr1", "chr2", "startBp1", "endBp1", "startBp2", "endBp2", "orient", "Species"),
  names_from = "rat_type",
  values_from = "value"
)

RATIOS_WIDE <-as.data.frame(RATIOS_WIDE)


NErat <- melted1[melted1$Var1=="Cbren",c(1,2)]
NErat$ratio_Ne <- 2
NErat[NErat$Var2=="Crema",]$ratio_Ne <-1
NErat[NErat$Var2=="Cinop",]$ratio_Ne <-1.3
NErat$ratio_Nestr <-0
NErat[NErat$Var2=="Crema",]$ratio_Nestr <- 1
NErat[NErat$Var2=="Cinop",]$ratio_Nestr <- 1
NErat[NErat$Var2=="Cnigo",]$ratio_Nestr <- 1
colnames(NErat) <-c("genome1","genome2","ratio_Ne","ratio_Ne_str")
RATIOS_WIDEP <-merge(RATIOS_WIDE,phylo, by =c("genome1","genome2"),all.x=TRUE)

RATIOS_WIDEP <-merge(RATIOS_WIDEP,NErat, by =c("genome1","genome2"),all.x=TRUE)





RATIOS_REP2 <- RATIOS_WIDEP %>%
  mutate(replication_factor = floor((endBp1 - startBp1) / 10000) + 1) %>%
  uncount(replication_factor, .remove = FALSE)


############## lets' see the effect size and signif
############## of differences between species using the replicated dataframe

library(dplyr)
library(broom)

RATIOS_REP2_only_ratios <-  RATIOS_REP2 %>% select(Species, starts_with("ratio"))
RATIOS_REP2_only_ratios <- RATIOS_REP2_only_ratios[,-c(9,10)]

# across(starts_with("ratio"), list(mean = mean, sd = sd, median=median, tstat=~t.test(.x, mu = 1)$statistic, pval= ~t.test(.x, mu = 1)$p.value))) %>% as.data.frame()


result_ratios_species <- RATIOS_REP2_only_ratios %>%
  group_by(Species) %>%
  summarize(
    across(starts_with("ratio"), list(mean = mean, sd = sd, median=median, wxstat=~wilcox.test(.x - 1, mu = 0, alternative = "two.sided")$statistic, pval=~wilcox.test(.x - 1, mu = 0, alternative = "two.sided")$p.value))) %>% as.data.frame()

#wilcox.test(x - 1, mu = 0, alternative = "two.sided")

write.table(result_ratios_species, file = "Table_syntenic_ratios_mean_sd_wxtest.txt", sep = "\t", quote=F, row.names = FALSE)
############models


#model Ne_ratio ~ phylo + gene_number + block_size + intron_size + exon_size + intergenic

RATIOS_REP3<- RATIOS_REP2[!(RATIOS_REP2$genome2 %in% c("Cnigo","Cbrig")),]
CnigoCbrig <- RATIOS_REP2[(RATIOS_REP2$genome2 %in% c("Cnigo","Cbrig")),]
set.seed(123)






formula3 <- ratio_Ne_str ~ ratio_exon_length + ratio_intron_length +
  ratio_uncovered_region + ratio_genecount + ratio_exoncount +
  ratio_introncount + ratio_size + ratio_size:ratio_genecount + ratio_exon_length:ratio_intron_length + ratio_exoncount:ratio_introncount + ratio_introncount:ratio_intron_length + ratio_exoncount:ratio_exon_length +

  formula4 <- ratio_Ne_str ~ ratio_exon_length + ratio_intron_length +
  ratio_uncovered_region + ratio_genecount + ratio_exoncount +
  ratio_introncount + ratio_size + ratio_size:ratio_genecount + ratio_exon_length:ratio_intron_length + ratio_exoncount:ratio_introncount + ratio_introncount:ratio_intron_length + ratio_exoncount:ratio_exon_length  + ratio_size:ratio_uncovered_region + ratio_uncovered_region:ratio_genecount




logistic_model3 <- glm(formula3, data = RATIOS_REP3, family = "binomial")
logistic_pred4 <- predict(logistic_model3, newdata = CnigoCbrig, type = "response")

logistic_accuracy4 <- sum((logistic_pred4 > 0.5) == CnigoCbrig$ratio_Ne_str) / nrow(CnigoCbrig)
logistic_accuracy4
#[1] 0.6791565

set.seed(123)

#nope
#normalized_data <- as.data.frame(scale(RATIOS_REP3[, c(11:18)]))
#outliers <- apply(abs(normalized_data) > 3, 1, any)
#cleaned_data <- normalized_data[!outliers, ]
#cleaned_data$ratio_Ne_str <- RATIOS_REP3[!outliers,20]

#logistic_model4 <- glm(formula4, data = cleaned_data, family = "binomial", maxit = 100)
#logistic_pred5 <- predict(logistic_model4, newdata = CnigoCbrig[, c(11:18,20)], type = "response")
#logistic_accuracy5 <- sum((logistic_pred5 > 0.5) == CnigoCbrig$ratio_Ne_str) / nrow(CnigoCbrig)





split_index <- createDataPartition(RATIOS_REP3$ratio_Ne, p = 0.7, list = FALSE)
train_set <- RATIOS_REP3[split_index, ]
#test_data <- RATIOS_REP3[-split_index, ]
valid_set <- RATIOS_REP3[-split_index, ]


dtrain <- xgb.DMatrix(data = as.matrix(train_set[, -c(1:10, 18:21)]), label = train_set$ratio_Ne_str)
dvalid <- xgb.DMatrix(data = as.matrix(valid_set[, -c(1:10, 18:21)]), label = valid_set$ratio_Ne_str)
dtest <- xgb.DMatrix(data = as.matrix(CnigoCbrig[, -c(1:10, 18:21)]), label = CnigoCbrig$ratio_Ne_str)

watchlist <- list(train = dtrain, valid = dvalid)


params <- list(
  objective = "binary:logistic",
  eta = 0.1,
  max_depth = 3,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  alpha = 0.01,
  lambda = 0.01
)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  watchlist = watchlist,
  eval_metric = "logloss",
  early_stopping_rounds = 10
)


#log <- xgb_model$evaluation_log

train_log <- xgb_model$evaluation_log$train_logloss
valid_log <- xgb_model$evaluation_log$valid_logloss

# Convert the logs to data frames
train_log_df <- as.data.frame(train_log)
valid_log_df <- as.data.frame(valid_log)
colnames(train_log_df) <- "logloss"
colnames(valid_log_df) <- "logloss"

# Combine the logs
combined_log_df <- rbind(
  transform(train_log_df, set = "train"),
  transform(valid_log_df, set = "valid")
)
combined_log_df$iter <- rep(c(1:100),2)
# Plot both training and validation loss over iterations
ggplot(combined_log_df, aes(x = iter, y = logloss, color = set)) +
  geom_line() +
  labs(title = "XGBoost Training and Validation Progress",
       x = "Iteration",
       y = "Log Loss",
       color = "Dataset") +
  theme_minimal()






xgb_pred3 <- predict(xgb_model, as.matrix(CnigoCbrig[, -c(1:10,18:21)]))
xgb_accuracy3 <- sum((xgb_pred3 > 0.5) == CnigoCbrig$ratio_Ne_str) / nrow(CnigoCbrig)
xgb_accuracy3
#[1] 0.6827525



#########
library(pROC)
library(caTools)
library(gridExtra)
library(gtable)

roc(CnigoCbrig$ratio_Ne_str, logistic_pred4,plot = TRUE, print.auc = TRUE)
#Area under the curve: 0.7363
roc(CnigoCbrig$ratio_Ne_str, xgb_pred3,plot = TRUE, print.auc = TRUE)
#Area under the curve: 0.6773

logistic_roc3 <- roc(CnigoCbrig$ratio_Ne_str, logistic_pred4,plot = TRUE, print.auc = TRUE)

prcoordglm<-coords(logistic_roc3, "all", ret = c("recall", "precision"))
prcoordglm <-prcoordglm[complete.cases(prcoordglm),]
trapz(prcoordglm$precision,prcoordglm$recall)
#[1] 0.2700225
xgb_roc3 <- roc(CnigoCbrig$ratio_Ne_str, xgb_pred3,plot = TRUE, print.auc = TRUE)
prcoordxgb<-coords(xgb_roc3, "all", ret = c("recall", "precision"))
prcoordxgb <-prcoordxgb[complete.cases(prcoordxgb),]
trapz(prcoordxgb$precision,prcoordxgb$recall)
#[1] 0.2275698


roc1 <- pROC::roc(CnigoCbrig$ratio_Ne_str, logistic_pred4, plot=FALSE,
                  legacy.axes=TRUE, percent=FALSE)

roc2 <- pROC::roc(CnigoCbrig$ratio_Ne_str, xgb_pred3,plot=FALSE,
                  legacy.axes=TRUE, percent=FALSE)


# compute recall and precision at each ROC curve's threshold
prcoords <- pROC::coords(roc1, "all", ret = c("threshold", "recall", "precision","sensitivity", "specificity"), transpose = FALSE)
prcoords2 <- pROC::coords(roc2, "all", ret = c("threshold", "recall", "precision","sensitivity", "specificity"), transpose = FALSE)


prcoords$Type <- "glm"
prcoords2$Type <- "xgboost"
prcoordscomb <- rbind(prcoords,prcoords2)


mae_logistic <- mean(abs(CnigoCbrig$ratio_Ne_str - logistic_pred4))
mae_xgb <- mean(abs(CnigoCbrig$ratio_Ne_str - xgb_pred3))

mae_logistic
#[1] 0.3545458
mae_xgb
#[1] 0.3866966


mse_logistic <- mean((CnigoCbrig$ratio_Ne_str - logistic_pred4)^2)
mse_xgb <- mean((CnigoCbrig$ratio_Ne_str - xgb_pred3)^2)

mse_logistic
#[1] 0.2323144
mse_xgb
#[1] 0.2700679

prcoordscomb <- prcoordscomb[complete.cases(prcoordscomb),]
prcoordscomb <- prcoordscomb[is.finite(prcoordscomb$threshold),]
prcoordscomb <- prcoordscomb[-510,]
f1_score_glm <- mean(2 * (prcoordscomb[prcoordscomb$Type=="glm",]$precision * prcoordscomb[prcoordscomb$Type=="glm",]$recall) / (prcoordscomb[prcoordscomb$Type=="glm",]$precision + prcoordscomb[prcoordscomb$Type=="glm",]$recall))
f1_score_xgb <- mean(na.omit(2 * (prcoordscomb[prcoordscomb$Type=="xgboost",]$precision * prcoordscomb[prcoordscomb$Type=="xgboost",]$recall) / (prcoordscomb[prcoordscomb$Type=="xgboost",]$precision + prcoordscomb[prcoordscomb$Type=="xgboost",]$recall)))

f1_score_glm
#[1] 0.5050161
f1_score_xgb
#[1] 0.511881



rocc <- ggplot(prcoordscomb, aes(1-sensitivity,specificity)) +
  geom_path(aes(color = Type)) +
  theme(aspect.ratio = 1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey80") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "grey40"),legend.key = element_blank(),plot.title = element_text(hjust = 0, face = "bold")) +
  ggtitle("B") + scale_color_manual(values =c("#1b9496", "#f0af18"),name="") +
  theme(legend.position = "top")

prc<- ggplot(prcoordscomb, aes(recall, precision)) +
  geom_path(aes(color = Type)) +
  theme(aspect.ratio = 1) +
  geom_abline(slope = 0, intercept = 0.5, linetype = "dotted", color = "grey80") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "grey40"),legend.key = element_blank(),plot.title = element_text(hjust = 0, face = "bold")) +
  ggtitle("A") + scale_color_manual(values =c("#1b9496", "#f0af18"),name="") +
  theme(legend.position = "top")



#library("knitr")
classStats <- data.frame(
  AUROC = c(0.7363, 0.6773),
  Accuracy = c(0.6792, 0.6827),
  MAE = c(0.3545, 0.3867),
  MSE = c(0.2323, 0.2701),
  F1 = c(0.505, 0.5119)
)


rownames(classStats) <- c("glm", "xgboost")
colnames(classStats) <- c("AUROC", "Accuracy", "MAE", "MSE", "F1")
classStats <- t(classStats)





tt1 <- ttheme_minimal(core=list(bg_params = list(fill=c(rep("#a5e4e6",5),rep("#ffe5ad",5)))))
table_grob2 <- tableGrob(classStats,theme= tt1)


title <- textGrob("C", gp = gpar(fontsize = 14, fontface = "bold"), hjust = 6, vjust = 0.5)
padding <- unit(50, "mm")

# add a row for the title
table2 <- gtable_add_rows(table_grob2, heights = grobHeight(title) + padding, pos = 0)
table2 <- gtable_add_grob(table2, title, 1, 1, 1, ncol(table2))
bottom_padding <- unit(61, "mm")
table2 <- gtable_add_rows(table2, heights = bottom_padding, pos = -1)



FIGS6 <- (prc + rocc +  plot_layout(guides = 'collect')  & theme(legend.position = 'bottom',plot.margin = margin(1,0.5,1,0.5) ) ) + table2

tiff(filename = "FigS6.tiff",width = 10.2,height = 6.2, res=300, units = "in") ### Figure S4
plot(FIGS6)
dev.off()

##########################

importance_values <- xgb.importance(model = xgb_model)
print(importance_values)
#Feature       Gain      Cover  Frequency Importance
#1: uncovered_region 0.40011176 0.31714116 0.26131387 0.40011176
#2:        genecount 0.16969083 0.17167893 0.15036496 0.16969083
#3:             size 0.15151022 0.15088232 0.15766423 0.15151022
#4:      exon_length 0.13152758 0.16468434 0.14890511 0.13152758
#5:    intron_length 0.07634120 0.09615725 0.14598540 0.07634120
#6:        exoncount 0.04341009 0.04293925 0.06131387 0.04341009
#7:      introncount 0.02740833 0.05651675 0.07445255 0.02740833

#importance_values$Feature <-gsub("ratio_","",importance_values$Feature)


coefs <- summary(logistic_model3)$coefficients
coefs <- coefs[order(-abs(coefs[, "Estimate"])), ]
#########
coefs <- summary(logistic_model3)$coefficients
coefs <- coefs[order(-abs(coefs[, "Estimate"])), ]

coefs <-as.data.frame(coefs)
#Estimate Std. Error     z value      Pr(>|z|)
#exoncount                  16.9561065 3.04740061   5.5641213  2.634767e-08
#introncount               -15.7180680 2.55025361  -6.1633353  7.122848e-10
#uncovered_region           14.1478671 0.30725704  46.0457057  0.000000e+00
#size                      -10.8729662 0.55634611 -19.5435287  4.682998e-85
#genecount                  -7.0760426 0.58665740 -12.0616268  1.684229e-33
#exon_length:exoncount       4.7776788 0.40125716  11.9067751  1.091163e-32
#exon_length                -4.6994373 0.49340307  -9.5245401  1.657744e-21
#intron_length               4.0598652 0.18991253  21.3775531 2.161704e-101
#intron_length:introncount  -2.6746286 0.16814764 -15.9064290  5.718351e-57
#(Intercept)                 1.6066372 0.06675662  24.0670839 5.530475e-128
#genecount:size             -1.1079709 0.11754264  -9.4261183  4.255462e-21
#exoncount:introncount       0.7664599 0.33139990   2.3127945  2.073394e-02
#exon_length:intron_length   0.1787865 0.18460249   0.9684946  3.327974e-01

coefs$variable <- rownames(coefs)

colnames(coefs) <-c("Estimate",   "Std.Error", "z_value",    "pval",   "variable")
coefs$padj <- p.adjust(coefs$pval, method = "bonferroni")



coefs$variable <- gsub("ratio_","",coefs$variable)
coefs$variable <-gsub("count"," count", coefs$variable)
coefs$variable <-gsub("_length"," length", coefs$variable)
coefs$variable <-gsub("uncovered_region", "intergenic",coefs$variable)
coefs$variable <-gsub("size","block size", coefs$variable)


ordered_levelsglm <- coefs$variable[order(coefs$z_value)]

coefs$variable <- factor(coefs$variable, levels = ordered_levelsglm)



coefplot <-ggplot(coefs, aes(x = variable, y = z_value)) +
  geom_bar(stat = "identity", fill = "#1b9496") +
  labs(title = "B",
       x = "",
       y = "Standardized Coefficient") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(hjust = 0, face = "bold"),
    legend.position = "top",
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text.y = element_text(face = "italic", angle = 0, hjust = 0),
    strip.text = element_text(size = 10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "grey40")
  ) +
  coord_flip() + geom_text(
    aes( label = ifelse(padj < 0.0001, "***",
                        ifelse(padj < 0.001, "**",
                               ifelse(padj < 0.01, "*",
                                      ifelse(padj < 0.05, ".", ""))))),
    vjust = 0.5,hjust=ifelse(coefs$z_value > 0,-0.05,1.05),
    size = 4,
    color = "grey25"
  )





###################################################



#install.packages("DALEX")
library(DALEX)

explainer <- explain(logistic_model3, x= CnigoCbrig$ratio_Ne_str, y = logistic_pred4, data = CnigoCbrig)

var_importanceglm <- variable_importance(explainer)

importance_dfGLM <- as.data.frame(var_importanceglm)
importance_dfGLM <-importance_dfGLM[importance_dfGLM$permutation==10,]

importance_dfGLM <- importance_dfGLM[grepl("ratio", importance_dfGLM$variable), ]
importance_dfGLM <- importance_dfGLM[!(grepl("Ne", importance_dfGLM$variable)), ]


importance_dfGLM$variable <- gsub("ratio_","",importance_dfGLM$variable)
importance_dfGLM$variable <-gsub("count"," count", importance_dfGLM$variable)
importance_dfGLM$variable <-gsub("_length"," length", importance_dfGLM$variable)
importance_dfGLM$variable <-gsub("uncovered_region", "intergenic",importance_dfGLM$variable)
importance_dfGLM$variable <-gsub("size","block size", importance_dfGLM$variable)


ordered_levels <- importance_dfGLM$variable[order(importance_dfGLM$dropout_loss)]

importance_dfGLM$variable <- factor(importance_dfGLM$variable, levels = ordered_levels)



imppl <- ggplot(importance_dfGLM, aes(x = variable, y = dropout_loss)) +
  geom_bar(stat = "identity", fill = "#1b9496") +
  labs(title = "C",
       x = "",
       y = "Permutation Feature Importance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(hjust = 0, face = "bold"),
    legend.position = "top",
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text.y = element_text(face = "italic", angle = 0, hjust = 0),
    strip.text = element_text(size = 10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "grey40")
  ) +
  coord_flip()






#levels(as.factor(RATIOS$Type))
#[1] "block size"    "gene count"    "exon count"    "intron count"  "exon length"   "intron length"
#[7] "intergenic"

RATIOS$Type <- factor(RATIOS$Type, levels=c("intergenic", "block size",    "gene count",    "exon count",    "intron count",  "exon length",   "intron length"))


ratioexinInter<- ggplot(RATIOS[RATIOS$Type %in% c("intergenic","exon count"),], aes(x = startBp1, y = value, ordered = FALSE, col=Type,fill=Type)) + facet_grid(Species ~ chr1, scales="free_x") + theme_bw(base_size = 12) + labs(title = "A", x = expression(~ italic("C. brenneri") ~ " genome position (Mb)"), y = expression("Ratio " ~ italic("C.sp") ~"/" ~italic("C. brenneri"))) + theme(
  axis.text.x = element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.title = element_text(face = "bold", hjust = 0)
) + theme(legend.position = "top") + scale_colour_manual(values = c("#808687","#de7702"), name="") + scale_fill_manual(values = c("#808687", "#de7702"),name="")  + theme(
  strip.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))
)  +  theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12)) + scale_x_continuous(limits = c(0,NA),
                                                                                                                      labels = function(x)
                                                                                                                        x / 1000000
) + coord_cartesian(clip = "off")  + geom_hline(yintercept = log2(1), linetype = "dotted", linewidth=0.9, color = "black") + geom_rect(aes(xmin = startBp1, xmax = endBp1, ymin=log2(value)-0.2, ymax=log2(value)+0.2), linewidth=1.35, alpha = 0.75, color='NA') + theme(strip.text.y = element_text(face = "italic")) + ylim(-4,4)





FIG7<- (ratioexinInter / (coefplot+imppl)) + plot_layout( heights = c(3,1))
tiff(filename = "Fig7.tiff",width = 12.2,height = 10.9, res=300, units = "in") ## Figure 5
plot(FIG7)
dev.off()




#################### probability of fixation ## Lynch 2002
#N=10^2
#s=0.001
#fixpr = 2*N*2*s/(exp(2*N*s) - 1)


calculate_fixpr <- function(N, s) {
  return (2 * N * 2 * s / (exp(2 * N * s) - 1))
}

# Ne and s combinations
#N_values <- 10^seq(2, 7, by = 1)
#s_values <- 10^seq(1, -9, by = -1)

#more deatails
N_values <- 10^(seq(2, 8, by = 0.25))
#s_values <- 10^(seq(1, -9, by = -0.25))
s_values <- 10^(seq(1, -9, by = -0.25))


resultsFIXPR <- data.frame(N = numeric(0), s = numeric(0), fixpr = numeric(0))


for (N in N_values) {
  for (s in s_values) {
    fixpr <- calculate_fixpr(N, s)
    resultsFIXPR <- rbind(resultsFIXPR, data.frame(N = N, s = s, Ns=N*s, fixpr = fixpr))
  }
}


custom_palette <- colorRampPalette(c("black","#4fa0b8","#f0b326"))(100)

#  scale_fill_viridis_c() +
vertical_lines <- data.frame(x_intercept = c(10^3, 10^6, 10^7),
                             label = c("C. tropicalis / C. elegans / C. briggsae",
                                       "C. remanei",
                                       "C. brenneri"))


FIXPRplo <- ggplot(resultsFIXPR, aes(x = N, y = s, fill = fixpr)) +
  geom_tile() +
  scale_x_log10(breaks = 10^(2:8), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = 10^(-9:1), labels = scales::trans_format("log10", scales::math_format(10^.x, format = "%g"))) +
  scale_fill_gradientn(colors = custom_palette) +
  labs(title = "",
       x = "Population size (Ne)",
       y = "Selection coefficient for introns (s)",
       fill = "fixing\nprobability") +
  theme_bw(base_size = 12)  + theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))
  )  +  theme(strip.text = element_text(size = 12))  + theme(strip.text = element_text(size = 12)) +
  geom_vline(data = vertical_lines, aes(xintercept = x_intercept), linetype = "dashed", color = "white",size=1.5) +
  geom_text(data = vertical_lines, aes(x = x_intercept, label = label), color = "white", angle = 90, vjust = 1.5, hjust = 0,size =4.5,  fontface = "italic")



tiff(filename = "FigS_prob.tiff",width = 7.2,height = 5.2, res=300, units = "in") ## Figure 1
plot(FIXPRplo)
dev.off()
