library(tidyverse)
library(maftools)
library(ggpubr)

setwd("/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/")
FGR <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/FGR_somatic.maf")
df_FGR <- FGR@data %>% 
  filter(dp > 20)

Normal <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/Normal_somatic.maf")
df_Normal <- Normal@data %>% 
  filter(dp > 20)

# Mutations---------------------------------------------------------------------
table(df_FGR$Tumor_Sample_Barcode)
table(df_Normal$Tumor_Sample_Barcode)

df_mutations <- data.frame(samples = c(names(table(df_FGR$Tumor_Sample_Barcode)), names(table(df_Normal$Tumor_Sample_Barcode))),
                 mutations = c(table(df_FGR$Tumor_Sample_Barcode), table(df_Normal$Tumor_Sample_Barcode)))
df_mutations <- df_mutations %>% 
  mutate(group = ifelse(grepl("N", samples), "Normal", "FGR"),
         twins = ifelse(grepl("F1", samples), "F1", "F2"))
df_mutations <- df_mutations %>% 
  group_by(twins, group) %>% 
  summarise(samples = samples, mutations = mutations, N = n(), mean = mean(mutations), sd = sd(mutations), se = sd(mutations)/sqrt(N))

ggplot(df_mutations, aes(x = group, y = mutations, fill = twins)) +
  geom_point(position = position_dodge(0.9), size = 1) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Group", y = "Mutation burden", title = "Somatic mutations number") +
  scale_fill_manual(values = c("#8093f1","#72ddf7")) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5))

ggplot(df_mutations, aes(x = group, y = mutations, color = group)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
               geom = "pointrange", color = "black", size = 0.5) +
  ylim(0, 400) +
  labs(x = "Group", y = "Mutation burden") +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "two.sided"), size = 5) +
  theme_classic()

# Median VAF--------------------------------------------------------------------
df_vaf <- rbind(df_FGR, df_Normal)
df_vaf <- df_vaf %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  summarise(median_vaf = median(vaf))
df_vaf <- df_vaf %>% 
  mutate(group = ifelse(grepl("N", Tumor_Sample_Barcode), "Normal", "FGR"),
         twins = ifelse(grepl("F1", Tumor_Sample_Barcode), "F1", "F2")) %>% 
  group_by(twins, group) %>% 
  summarise(samples = Tumor_Sample_Barcode, median_vaf = median_vaf, N = n(), mean = mean(median_vaf), sd = sd(median_vaf), se = sd(median_vaf)/sqrt(N))


ggplot(df_vaf, aes(x = group, y = median_vaf, fill = twins)) +
  geom_point(position = position_dodge(0.9), size = 1) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Group", y = "Median VAF", title = "VAF") +
  scale_fill_manual(values = c("#8093f1","#72ddf7")) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5))

ggplot(df_vaf, aes(x = group, y = median_vaf, color = group)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
               geom = "pointrange", color = "black", size = 0.5) +
  ylim(0, 0.6) +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = "two.sided"), size = 5) +
  labs(x = "Group", y = "Median VAF") +
  theme_classic()

# vaf distribution -------------------------------------------------------------
df_all <- rbind(df_FGR, df_Normal)
df_all <- df_all %>% 
  mutate(vaf_bin = cut(vaf, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE, labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")))
df_summary <- df_all %>%
  group_by(Tumor_Sample_Barcode, vaf_bin) %>%
  summarise(count = n())
ggplot(df_summary, aes(x = vaf_bin, y = count, fill = vaf_bin)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Tumor_Sample_Barcode) +
  theme_classic()

# intersection -----------------------------------------------------------------
df_FGR_venn <- df_FGR %>% 
  mutate(var = paste(Chromosome, Start_Position, Strand, Reference_Allele, Tumor_Seq_Allele2, sep = "_"),
         sample = str_split(Tumor_Sample_Barcode, "-", simplify = T)[, 2],
         twins = str_split(Tumor_Sample_Barcode, "-", simplify = T)[, 1])
df_Normal_venn <- df_Normal %>% 
  mutate(var = paste(Chromosome, Start_Position, Strand, Reference_Allele, Tumor_Seq_Allele2, sep = "_"),
         sample = str_split(Tumor_Sample_Barcode, "_", simplify = T)[, 2],
         twins = str_split(Tumor_Sample_Barcode, "_", simplify = T)[, 3])
df_FGR_venn <- df_FGR_venn %>% 
  group_by(sample, twins) %>% 
  summarise(vars = list(unique(var))) %>%
  pivot_wider(names_from = twins, values_from = vars)
df_Normal_venn <- df_Normal_venn %>% 
  group_by(sample, twins) %>% 
  summarise(vars = list(unique(var))) %>%
  pivot_wider(names_from = twins, values_from = vars)
df_FGR_venn <- df_FGR_venn %>%
  mutate(intersection_count = map2_int(F1, F2, ~ length(intersect(.x, .y))))
df_Normal_venn <- df_Normal_venn %>%
  mutate(intersection_count = map2_int(F1, F2, ~ length(intersect(.x, .y))))

library(ggVennDiagram)
library(gridExtra)
# FGR
plots <- list()
for (i in 1:nrow(df_FGR_venn)){
  sample <- df_FGR_venn$sample[i]
  data <- list(F1 = df_FGR_venn$F1[[i]], F2 = df_FGR_venn$F2[[i]])
  plots[[i]] <- ggVennDiagram(data) +
    coord_flip() +
    ggtitle(sample) +
    scale_fill_gradient(low="white",high = "#b9292b",name = "count") +
    theme(plot.title = element_text(size = 15,hjust = 0.5))
}
grid.arrange(grobs=plots, ncol=4, nrow=4)
# Normal
plots <- list()
for (i in 1:nrow(df_Normal_venn)){
  sample <- df_Normal_venn$sample[i]
  data <- list(F1 = df_Normal_venn$F1[[i]], F2 = df_Normal_venn$F2[[i]])
  plots[[i]] <- ggVennDiagram(data) +
    coord_flip() +
    ggtitle(sample) +
    scale_fill_gradient(low="white",high = "#b9292b",name = "count") +
    theme(plot.title = element_text(size = 15,hjust = 0.5))
}
grid.arrange(grobs=plots, ncol=3, nrow=2)

# vaf 0-0.2 distribution -------------------------------------------------------
df_all <- rbind(df_FGR, df_Normal)
df_all_vaf_020 <- df_all %>% 
  filter(vaf<=0.20)
df_all_vaf_020 <- df_all_vaf_020 %>% 
  mutate(vaf_bin = cut(vaf, breaks = c(0, 0.01, 0.02, 0.04, 0.10, 0.20), labels = c("0-0.01", "0.01-0.02", "0.02-0.04", "0.04-0.10", "0.10-0.20")))
df_summary <- df_all_vaf_020 %>%
  group_by(Tumor_Sample_Barcode, vaf_bin) %>%
  summarise(count = n())
ggplot(df_summary, aes(x = vaf_bin, y = count, fill = vaf_bin)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Tumor_Sample_Barcode) +
  theme_classic()

# intersection(vaf < 0.2) ------------------------------------------------------
df_FGR_vaf_020 <- df_FGR %>% 
  filter(vaf<=0.20)
df_Normal_vaf_020 <- df_Normal %>% 
  filter(vaf<=0.20)
df_FGR_venn <- df_FGR_vaf_020 %>% 
  mutate(var = paste(Chromosome, Start_Position, Strand, Reference_Allele, Tumor_Seq_Allele2, sep = "_"),
         sample = str_split(Tumor_Sample_Barcode, "-", simplify = T)[, 2],
         twins = str_split(Tumor_Sample_Barcode, "-", simplify = T)[, 1])
df_Normal_venn <- df_Normal_vaf_020 %>% 
  mutate(var = paste(Chromosome, Start_Position, Strand, Reference_Allele, Tumor_Seq_Allele2, sep = "_"),
         sample = str_split(Tumor_Sample_Barcode, "_", simplify = T)[, 2],
         twins = str_split(Tumor_Sample_Barcode, "_", simplify = T)[, 3])
df_FGR_venn <- df_FGR_venn %>% 
  group_by(sample, twins) %>% 
  summarise(vars = list(unique(var))) %>%
  pivot_wider(names_from = twins, values_from = vars)
df_Normal_venn <- df_Normal_venn %>% 
  group_by(sample, twins) %>% 
  summarise(vars = list(unique(var))) %>%
  pivot_wider(names_from = twins, values_from = vars)
df_FGR_venn <- df_FGR_venn %>%
  mutate(intersection_count = map2_int(F1, F2, ~ length(intersect(.x, .y))))
df_Normal_venn <- df_Normal_venn %>%
  mutate(intersection_count = map2_int(F1, F2, ~ length(intersect(.x, .y))))

library(ggVennDiagram)
library(gridExtra)
# FGR
plots <- list()
for (i in 1:nrow(df_FGR_venn)){
  sample <- df_FGR_venn$sample[i]
  data <- list(F1 = df_FGR_venn$F1[[i]], F2 = df_FGR_venn$F2[[i]])
  plots[[i]] <- ggVennDiagram(data) +
    coord_flip() +
    ggtitle(sample) +
    scale_fill_gradient(low="white",high = "#b9292b",name = "count") +
    theme(plot.title = element_text(size = 15,hjust = 0.5))
}
grid.arrange(grobs=plots, ncol=4, nrow=4)
# Normal
plots <- list()
for (i in 1:nrow(df_Normal_venn)){
  sample <- df_Normal_venn$sample[i]
  data <- list(F1 = df_Normal_venn$F1[[i]], F2 = df_Normal_venn$F2[[i]])
  plots[[i]] <- ggVennDiagram(data) +
    coord_flip() +
    ggtitle(sample) +
    scale_fill_gradient(low="white",high = "#b9292b",name = "count") +
    theme(plot.title = element_text(size = 15,hjust = 0.5))
}
grid.arrange(grobs=plots, ncol=3, nrow=2)

# SBS --------------------------------------------------------------------------
df_all_vaf_001 <- df_all %>% 
  filter(vaf<0.01)
df_all_vaf_002 <- df_all %>% 
  filter(0.01<=vaf & vaf<0.02)
df_all_vaf_004 <- df_all %>% 
  filter(0.02<=vaf & vaf<0.04)
df_all_vaf_020 <- df_all %>% 
  filter(0.04<=vaf & vaf<0.20)
write.table(df_all, "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/all_somatic_dp20.maf", quote = F, row.names = F, sep = "\t")
write.table(df_all_vaf_001, "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/all_somatic_dp20_vaf001.maf", quote = F, row.names = F, sep = "\t")
write.table(df_all_vaf_002, "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/all_somatic_dp20_vaf002.maf", quote = F, row.names = F, sep = "\t")
write.table(df_all_vaf_004, "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/all_somatic_dp20_vaf004.maf", quote = F, row.names = F, sep = "\t")
write.table(df_all_vaf_020, "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/all_somatic_dp20_vaf020.maf", quote = F, row.names = F, sep = "\t")
