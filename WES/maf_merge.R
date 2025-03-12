library(tidyverse)
library(maftools)


setwd("/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/")
FGR_maf_files <- list.files(path = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/", pattern = "^F.*vep.maf$", full.names = T)
FGR_maf_files
Normal_maf_files <- list.files(path = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/", pattern = "^N.*vep.maf$", full.names = T)
Normal_maf_files
FGR_mymaf <- maftools::merge_mafs(mafs = FGR_maf_files)
maftools::write.mafSummary(maf = FGR_mymaf, basename = "FGR_placentas")
Normal_mymaf <- maftools::merge_mafs(mafs = Normal_maf_files)
maftools::write.mafSummary(maf = Normal_mymaf, basename = "Normal_placentas")

FGR_maf <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/FGR_placentas_maftools.maf")
df_FGR <- FGR_maf@data
table(df_FGR$Source_MAF)

df_FGR_filter <- df_FGR %>% 
  mutate(vaf = t_alt_count/(t_alt_count + t_ref_count)) %>%
  mutate(dp = (t_alt_count + t_ref_count)) %>% 
  filter(FILTER == "PASS")
table(df_FGR_filter$Source_MAF)
write.table(df_FGR_filter, "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/FGR_somatic.maf", quote = F, row.names = F, sep = "\t")

# FGR <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/FGR_somatic.maf")
# FGR@variants.per.sample
# oncoplot(maf = FGR, top = 10)
# FGR.titv <- titv(maf = FGR, plot = F, useSyn = T)
# plotTiTv(res = FGR.titv)
# rainfallPlot(maf = FGR, detectChangePoints = T, pointSize = 0.4)
# laml.mutload = tcgaCompare(maf = FGR, cohortName = 'placenta_FGR', logscale = TRUE, capture_size = 20)
# plotVaf(maf = FGR)
# FGR.sig <- oncodrive(maf = FGR, minMut = 5, pvalMethod = 'zscore')
# plotOncodrive(res = FGR.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 2)
# 
# library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
# FGR.tnm = trinucleotideMatrix(maf = FGR, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
# FGR.tnm$nmf_matrix <- FGR.tnm$nmf_matrix + 0.001
# library(NMF)
# laml.sign = estimateSignatures(mat = FGR.tnm, nTry = 10, pConstant = 0.1)
# plotCophenetic(res = laml.sign)
# FGR.sig <- extractSignatures(mat = FGR.tnm, n = 3)
# FGR.og30.cosm = compareSignatures(nmfRes = FGR.sig, sig_db = "legacy")
# FGR.v3.cosm = compareSignatures(nmfRes = FGR.sig, sig_db = "SBS")
# FGR.v3.cosm
# maftools::plotSignatures(nmfRes = FGR.sig, title_size = 1.2, sig_db = "SBS")


Normal_maf <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/Normal_placentas_maftools.maf")
df_Normal <- Normal_maf@data
df_Normal_filter <- df_Normal %>% 
  mutate(vaf = t_alt_count/(t_alt_count + t_ref_count)) %>%
  mutate(dp = (t_alt_count + t_ref_count)) %>% 
  filter(FILTER == "PASS")
table(df_Normal_filter$Source_MAF)
write.table(df_Normal_filter, "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/Normal_somatic.maf", quote = F, row.names = F, sep = "\t")

# Normal <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/Normal_somatic.maf")
# Normal@variants.per.sample
# oncoplot(maf = Normal, top = 10)
# Normal.titv <- titv(maf = Normal, plot = F, useSyn = T)
# plotTiTv(res = Normal.titv)
# rainfallPlot(maf = Normal, detectChangePoints = T, pointSize = 0.4)
# laml.mutload = tcgaCompare(maf = Normal, cohortName = 'placenta_FGR', logscale = TRUE, capture_size = 20)
# plotVaf(maf = Normal)
# 
# library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
# Normal.tnm = trinucleotideMatrix(maf = Normal, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
# Normal.tnm
# Normal.sig <- extractSignatures(mat = Normal.tnm, n = 3)
# maftools::plotSignatures(nmfRes = Normal.sig, title_size = 1.2, sig_db = "SBS")
# 
# mafCompare(m1 = FGR, m2 = Normal)

df_all_unfilter <- rbind(df_FGR, df_Normal)
df_all_unfilter <- df_all_unfilter %>% 
  mutate(vaf = t_alt_count/(t_alt_count + t_ref_count)) %>%
  mutate(dp = (t_alt_count + t_ref_count)) %>% 
  filter(dp > 20)

df_all_unfilter_C2A <- df_all_unfilter %>% 
  filter(Reference_Allele == "C" & Tumor_Seq_Allele2 == "A")
df_all_unfilter_C2T <- df_all_unfilter %>% 
  filter(Reference_Allele == "C" & Tumor_Seq_Allele2 == "T")
sort(table(df_all_unfilter_C2A$FILTER))
sort(table(df_all_unfilter_C2T$FILTER))


