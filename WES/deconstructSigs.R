library(tidyverse)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)

FGR <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/FGR_somatic.maf")
df_FGR <- FGR@data %>% 
  filter(dp > 20) %>% 
  filter(vaf > 0.05 & vaf < 0.95)

Normal <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/Normal_somatic.maf")
df_Normal <- Normal@data %>% 
  filter(dp > 20) %>% 
  filter(vaf > 0.05 & vaf < 0.95)

chrom <- c(paste0("chr",1:22), "chrX", "chrY", "chrM")
df_FGR <- filter(df_FGR, Chromosome %in% chrom)
df_Normal <- filter(df_Normal, Chromosome %in% chrom)
unique(df_FGR$Chromosome)
unique(df_Normal$Chromosome)

# FGR --------------------------------------------------------------------------
sigs.input <- mut.to.sigs.input(mut.ref = df_FGR, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

sigs.input

# Determine the signatures contributing an already normalized sample
par(mfrow=c(5, 6),mar=c(2, 2, 2, 2))
for(sample in unique(df_FGR$Tumor_Sample_Barcode)){
  test = whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample,
                         contexts.needed = TRUE)
  makePie(test)
}
dev.off()

# Normal --------------------------------------------------------------------------
sigs.input <- mut.to.sigs.input(mut.ref = df_Normal, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

sigs.input

# Determine the signatures contributing an already normalized sample
par(mfrow=c(3, 4),mar=c(2, 2, 2, 2))
for(sample in unique(df_Normal$Tumor_Sample_Barcode)){
  test = whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample,
                         contexts.needed = TRUE)
  makePie(test)
}
dev.off()