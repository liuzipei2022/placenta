library(tidyverse)

# 创建数据框
df <- read.table("/picb/lilab5/liuzipei/placenta/WES/3.bam_qc/all_samples_summary.txt", sep = "\t", header = T)

df <- df %>% 
  select(Sample, Depth) %>% 
  mutate(twins = ifelse(grepl("F1", Sample), "F1", "F2"),
         group = ifelse(grepl("N", Sample), "Normal", "FGR")) %>% 
  group_by(twins, group) %>%
  summarise(Sample = Sample, Depth = Depth, N = n(), mean = mean(Depth), sd = sd(Depth), se = sd(Depth)/sqrt(N))
  

ggplot(df, aes(x = group, y = Depth, fill = twins)) +
  geom_point(position = position_dodge(0.9), size = 1) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Group", y = "Average depth", title = "WES") +
  scale_fill_manual(values = c("#8093f1","#72ddf7")) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5)) 


