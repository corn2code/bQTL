#this code will input the output from LD analysis

library(data.table)

LD.R <- fread("bQTL.cisEQTL.10kb.LD.Random", data.table = F)
colnames(LD.R) <- c("CHR_A", "BP_A", "SNP_A",  "CHR_B", "BP_B", "SNP_B", "R2")

#remove comparison in both initial datasets
LD.R.1 <- LD.R[!(LD.R$SNP_A %in% bQTL.random[,4] & LD.R$SNP_B %in% bQTL.random[,4]) & # SNPs used in bQTL calculation in eGWAS340
                         !(LD.R$SNP_A %in% Sun.cis[,15] & LD.R$SNP_B %in% Sun.cis[,15]),] # Sun.cis = cis-eQTL SNPs

#select comparisons with R2 greater than 0.6
LD.R.1.6 <- LD.R.1[which(LD.R.1$R2 > 0.6),]

LD.true <- fread("bQTL.cisEQTL.10kb.LD", data.table = F)
colnames(LD.true) <- c("CHR_A", "BP_A", "SNP_A",  "CHR_B", "BP_B", "SNP_B", "R2")

#remove comparison in both initial datasets
LD.true.1 <- LD.true[!(LD.true$SNP_A %in% bQTL[,1] & LD.true$SNP_B %in% bQTL[,1]) & #bQTL = bQTL in eGWAS340
                               !(LD.true$SNP_A %in% Sun.cis[,15] & LD.true$SNP_B %in% Sun.cis[,15]),] # Sun.cis = cis-eQTL SNPs

#select comparisons with R2 greater than 0.6
LD.true.1.6 <- LD.true.1[which(LD.true.1$R2 > 0.6),]


write.csv(LD.true.1.6,"LD.bQTL.cisEQTL.csv", row.names = F)
write.csv(LD.R.1.6,"LD.ControlBQTL.cisEQTL.csv", row.names = F)