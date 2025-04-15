setwd("/work/schnablelab/vladimir/eQTL_bQTL/Revision1")


dist.10kb <- fread("bQTL_WW_SNP.clumMetGeno.GenoOnly_v2.version4.genes.txt", data.table = F)
dist.10kb$V4.num <- as.numeric(dist.10kb$V4)
#dist.10kb$distance <-  dist.10kb$V4
dist.10kb$V4.num[which(dist.10kb$V4 == "intergenic")] <- 11000
dist.10kb$V4.num[which(dist.10kb$V4 == "intragenic")] <- -1000

#add col names
colnames(dist.10kb) <- c("bQTL","start","end","distance","target","dist.num")


eGWAS_snps.txt <- fread("eGWAS_snps.txt", data.table = F)
colnames(eGWAS_snps.txt) <- c("rs", "alleles", "chr","pos")

dist.10kb.eGWAS <- merge(dist.10kb, eGWAS_snps.txt, by.x = "end", by.y = "pos")


# load suns eQTLs, rs and pos are different
Sun <- fread("e340_eSummary_filt.csv", data.table = F)
Sun <- Sun[which(Sun$Type == "Cis"),]
NROW(unique(Sun$gene))
NROW(unique(Sun$snps))

#
length(intersect(Sun$snps, dist.10kb.eGWAS$rs)) # 94 similar 

SNPS.bQTL <- data.frame(c(Sun$snps, dist.10kb.eGWAS$rs))

#save list of snps
write.table(SNPS.bQTL[,1],paste0("LD/SNPS.bQTL.txt"),row.names=FALSE,sep="\t", quote = FALSE, col.names = F)

system(paste0("wc -l ",paste0("LD/SNPS.bQTL.txt")))

# Construct the PLINK command you want to run as a string

plinkPath <- "/util/opt/anaconda/deployed-conda-envs/packages/plink/envs/plink-1.90b4/bin/plink"
plink_command <- paste0(plinkPath," --bfile eGWAS340_maf05het02 --r2 --ld-window-r2 0 --ld-window-kb 10 --ld-window 99999  --extract LD/SNPS.bQTL.txt  --out LD/bQTL.cisEQTL.10kb.LD")

#plink --bfile eGWAS340_maf05het02 --r2 --ld-window-r2 0 --ld-window-kb 10 
#--ld-window 99999  --extract LD/SNPS.bQTL.txt  --out LD/bQTL.cisEQTL.10kb.LD


system(plink_command)


#load plink results
LD.true <- fread(paste0("LD/bQTL.cisEQTL.10kb.LD.ld"),data.table = F)
colnames(LD.true) <- c("CHR_A", "BP_A", "SNP_A",  "CHR_B", "BP_B", "SNP_B", "R2")

#remove comparison in both initial datasets
LD.true.1 <- LD.true[!(LD.true$SNP_A %in% dist.10kb.eGWAS[,7] & LD.true$SNP_B %in% dist.10kb.eGWAS[,7]) & #bQTL = bQTL in eGWAS340
                       !(LD.true$SNP_A %in% Sun[,15] & LD.true$SNP_B %in% Sun[,15]),] # Sun.cis = cis-eQTL SNPs

#select comparisons with R2 greater than 0.6
LD.true.1.6 <- LD.true.1[which(LD.true.1$R2 > 0.6),] # 1775
nrow(LD.true.1.6)

fwrite(LD.true.1.6,"LD.bQTL.cisEQTL.csv")

###### now put together a random dataset to use as a control
## first look at the distribution 
dist.10kb

distance <-  c("intragenic","0-1", "1-2", "2-3","3-4","4-5", "5-6", "6-7","7-8",
               "8-9", "9-10", "intergenic") # h$mids coming from hist 
#mpg moves the lables
svg("hist_distance.38291.bQTL.svg", height = 6, width = 8)
#### find the distribution of the distance
par(mar=c(4,4,4,4)+0.1,mgp=c(3,1,0)) # sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
h <- hist(dist.10kb.eGWAS$dist.num, breaks = 12, main = "bQTL.V4 in 12M SNP dataset used for cis-EQTL", xaxt = "n", xlab="distance to gene in KB")
axis(1, at=h$mids, labels = FALSE)
text(x = h$mids,
     ## Move labels to just below bottom of chart.
     y = -1900,
     ## Use names from the data list.
     labels = distance,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 1.2)
text(h$mids[1:11],h$counts[1:11],labels=h$counts[1:11], adj=c(0.1, -0.3), srt=45)
text(h$mids[12],h$counts[12],labels=h$counts[12], adj=c(1.5, 1), srt=0)
dev.off()

h$counts

#### load random data
M2 <- fread("sorted_random_2M_lines_v2.version4.genes.txt")
M2$V4.num <- as.numeric(M2$V4)
M2$V4.num[which(M2$V4 == "intergenic")] <- 11000
M2$V4.num[which(M2$V4 == "intragenic")] <- -1000
colnames(M2) <- c("chr","start","end","distance","target","dist.num")


M2.eGWAS <- merge(M2, eGWAS_snps.txt, by.x = "end", by.y = "pos")

bQTL.random <- data.frame()
set.seed(4)
bQTL.random <- M2.eGWAS[sample(which(M2.eGWAS$dist.num == -1000), 6481),] #!!! number are corrected as for may 11 2023
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 0 &
                                                                 M2.eGWAS$dist.num <= 1000), 8005),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 1001 &
                                                                 M2.eGWAS$dist.num <= 2000), 2846),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 2001 &
                                                                 M2.eGWAS$dist.num <= 3000), 1741),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 3001 &
                                                                 M2.eGWAS$dist.num <= 4000), 1135),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 4001 &
                                                                 M2.eGWAS$dist.num <= 5000), 914),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 5001 &
                                                                 M2.eGWAS$dist.num <= 6000), 723),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 6001 &
                                                                 M2.eGWAS$dist.num <= 7000), 621),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 7001 &
                                                                 M2.eGWAS$dist.num <= 8000), 460),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 8001 &
                                                                 M2.eGWAS$dist.num <= 9000), 435),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num >= 9001 &
                                                                 M2.eGWAS$dist.num <= 10000), 391),])
set.seed(4)
bQTL.random <- rbind(bQTL.random, M2.eGWAS[sample(which(M2.eGWAS$dist.num == 11000), 14539),])
NROW(bQTL.random)
bQTL.random <- as.data.frame(bQTL.random)

length(intersect(Sun$snps, bQTL.random$rs)) # 55 similar versus 94 similar with true 

SNPS.random <- data.frame(c(Sun$snps, bQTL.random$rs))

#save list of snps
write.table(SNPS.random[,1],paste0("LD/SNPS.random.txt"),row.names=FALSE,sep="\t", quote = FALSE, col.names = F)

system(paste0("wc -l ",paste0("LD/SNPS.random.txt")))

# Construct the PLINK command you want to run as a string

plinkPath <- "/util/opt/anaconda/deployed-conda-envs/packages/plink/envs/plink-1.90b4/bin/plink"
plink_command2 <- paste0(plinkPath," --bfile eGWAS340_maf05het02 --r2 --ld-window-r2 0 --ld-window-kb 10 --ld-window 99999  --extract LD/SNPS.random.txt  --out LD/random.cisEQTL.10kb.LD")

system(plink_command2)


#load plink results
LD.random<- fread(paste0("LD/random.cisEQTL.10kb.LD.ld"),data.table = F)
colnames(LD.random) <- c("CHR_A", "BP_A", "SNP_A",  "CHR_B", "BP_B", "SNP_B", "R2")

#remove comparison in both initial datasets
LD.random.1 <- LD.random[!(LD.random$SNP_A %in% bQTL.random[,7] & LD.random$SNP_B %in% bQTL.random[,7]) & #bQTL = bQTL in eGWAS340
                       !(LD.random$SNP_A %in% Sun[,15] & LD.random$SNP_B %in% Sun[,15]),] # Sun.cis = cis-eQTL SNPs

#select comparisons with R2 greater than 0.6
LD.random.1.6 <- LD.random.1[which(LD.random.1$R2 > 0.6),] # 1775
nrow(LD.random.1.6)

fwrite(LD.random.1.6,"LD.random.cisEQTL.csv")

# true / random
nrow(LD.true.1.6)/nrow(LD.random.1.6) # 1.28

nrow(LD.true.1.6)
nrow(LD.random.1.6)

1649/1232 # from before 1.338474

3478/2894


library(data.table)
library(ggplot2)
library(svglite)

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5),legend.title=element_text(size=10), 
             legend.text=element_text(size=9))


LD.true.1.6 <- fread("LD.bQTL.cisEQTL.csv",data.table=F)
LD.random.1.6 <- fread("LD.random.cisEQTL.csv",data.table=F)

intersect(LD.true.1.6$SNP_A, LD.true.1.6$SNP_B)
intersect(LD.random.1.6$SNP_A, LD.random.1.6$SNP_B)

# Combine the two datasets into one data frame with a new column to indicate the group
LD.true.1.6$group <- "bQTL "
LD.random.1.6$group <- "bgSNP"

combined_data <- rbind(LD.true.1.6, LD.random.1.6)

# Create the plot
ggplot(combined_data, aes(x = R2, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 9) +
  scale_fill_manual(values = c("blue", "grey")) +
  labs(x = "Linkage Disequilibrium ", y = "Frequency") +
  coord_flip()

ggsave("comparison_WW_SNP.clumMetGeno.GenoOnly.svg", height = 6, width = 6)
ggsave("comparison_WW_SNP.clumMetGeno.GenoOnly.png", height = 6, width = 6)

######

dta <- fread("LD.bQTL.cisEQTL.csv", data.table = F)
dta.R <- fread("LD.random.cisEQTL.csv", data.table = F)

bQTL <- hist(dta$R2)
bgSNP <- hist(dta.R$R2)

svg("LDcomparison.svg", height = 5, width = 7)
plot(bQTL, col="blue", xlim=c(.6,1), main = "", xlab = "Linkage Disequilibrium")  # first histogram
plot(bgSNP, col="grey", xlim=c(.6,1), add=T)  # second
legend("topright", legend=c("bQTL","bgSNP"), fill=c("blue","grey"))
dev.off()

