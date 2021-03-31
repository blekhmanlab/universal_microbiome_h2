############################################################################
#### Script 4: Humans and baboons share heritable taxa ####
############################################################################

library(ggplot2)
library(RColorBrewer)
library(lmerTest)
library(dplyr)

#### Comparing heritability estimates in humans vs baboons ####

#read in current baboon study, limit to only heritable traits
babs <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv", row.names = 1)
babs <- babs %>%
  mutate(baboon_h2_estimate = heritability) %>%
  mutate(baboon_is_heritable = ifelse(p_adjust < 0.1, "heritable", "NS")) %>%
  mutate(simple_name = phenotype) %>%
  filter(baboon_is_heritable == "heritable") %>%
  dplyr::select(baboon_h2_estimate, baboon_is_heritable, simple_name)

#read in processed human studies, limit to only heritable traits
humans <- read.csv("git_human_heritability_studies_noASV_or_unknown.csv", row.names = 1)
humans <- humans %>%
  filter(is_heritable == "heritable")

#limit datasets to only traits with exact name nmatches 
mods_graph <- merge(babs, humans, by = "simple_name")
mods_graph <- mods_graph %>%
  arrange(simple_name)

#add a descriptive label
mods_stats <- mods_graph %>%
  group_by(dataset) %>%
  dplyr::summarise(n_traits = n()) %>%
  mutate(`Human dataset and n shared heritable traits` = paste(dataset, " (", n_traits, ")", sep= ''))

mods_graph <- merge(mods_graph, mods_stats, by = "dataset")

# Linear mixed model, controlling for some taxa appearing in multiple datasets
lm1 <- lmer(h2_estimate ~ baboon_h2_estimate + (1|dataset), data = mods_graph)
summary(lm1)

# Calculate the correlation between mean human and baboon heritability
mean_h2 <- mods_graph %>%
  group_by(simple_name) %>%
  summarise(mean_human_h2 = mean(h2_estimate))
baboon_h2 <- mods_graph %>%
  dplyr::select(simple_name, baboon_h2_estimate) %>%
  distinct()
mean_h2 <- merge(mean_h2, baboon_h2, by = "simple_name")
cor.test(mean_h2$baboon_h2_estimate, mean_h2$mean_human_h2)


#### Fig. 2D ####

#specify colors
palEcol <-c('#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6')

t1<-theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),                               
          plot.background = element_blank(), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(size=.4), 
          axis.line.y = element_line(size=.4), 
          legend.position= "top",
          legend.text = element_text(color="black", size=12),
          legend.title = element_text(color="black", size=16),
          axis.title = element_text(color="black", size=20),
          axis.text = element_text(color="black", size=16)
)

Fig2D <- ggplot(mods_graph, aes(x=baboon_h2_estimate, y = h2_estimate)) + 
  geom_point(aes(color = `Human dataset and n shared heritable traits`), size = 3) +
  geom_smooth(method = 'lm', color = "black") +
  scale_color_manual(values = palEcol, guide = guide_legend(ncol = 2, title.position="top", title.hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.025, 0.05, 0.075, 0.1, 0.125), labels = c(0.025, 0.05, 0.075, 0.1, 0.125)) +
  labs(x = expression(paste("Baboon heritability (", h^2, ")"))) +
  labs(y = expression(paste("Human heritability (", h^2, ")"))) +
  t1