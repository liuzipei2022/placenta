library(tidyverse)
library(maftools)
library(ggpubr)

setwd("/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/")
FGR <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/FGR_somatic.maf")
df_FGR <- FGR@data %>% 
  filter(dp > 20) %>% 
  filter(vaf > 0.05 & vaf < 0.95)

Normal <- read.maf(maf = "/picb/lilab5/liuzipei/placenta/WES/5.somatic/filter/Normal_somatic.maf")
df_Normal <- Normal@data %>% 
  filter(dp > 20) %>% 
  filter(vaf > 0.05 & vaf < 0.95)

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

