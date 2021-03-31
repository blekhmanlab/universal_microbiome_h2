############################################################################
#### Script 6: Longitudinal sampling affects heritability estimation ####
############################################################################

library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(RColorBrewer)

#### Fig. 4A ####

m1 <- read.csv("git_metadata_and_community_phenotypes.csv", header = TRUE, row.names = 1)

#get it to like the highest 100 samples = 18 baboons
m2 <- m1 %>%
  dplyr::group_by(baboon_id) %>%
  dplyr::summarise(n_samples = n()) %>%
  filter(n_samples > 100) %>%
  droplevels()

m2 <- left_join(m2, m1, by = "baboon_id")

#add on tax
tax_all <- readRDS("git_CLR_table.RDS")
tax_all <- tax_all %>%
  dplyr::select(sample, f_Christensenellaceae)

m2 <- merge(m2, tax_all, by = "sample")

#order by who has the oldest sample
age_count <- m2 %>%
  dplyr::group_by(baboon_id) %>%
  dplyr::summarise(max_age = max(age)) %>%
  arrange(max_age)

age_count$age_order <- 1:nrow(age_count)

#order by modsel ordeer
age_count$baboon_id2 <- factor(age_count$baboon_id, levels = age_count$baboon_id[order(age_count$age_order)])
m2 <- merge(m2, age_count, by = "baboon_id")

#make Christensenellaceae numeric
m2$Christensenellaceae2 <- as.character(m2$f_Christensenellaceae)
m2$Christensenellaceae <- as.numeric(m2$Christensenellaceae2)

#theme and plot
t3<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major.x = element_blank(), 
  panel.grid.major.y = element_line(size=.4, colour = "grey", linetype = "dashed"), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line.x = element_line(size=.4), 
  axis.line.y = element_line(size=.4), 
  axis.ticks.y = element_blank(),
  legend.position = "top",
  axis.title = element_text(size = 20),
  axis.text.y=element_blank(),
  axis.text.x=element_text(size = 18),
  strip.text.y = element_blank(), #don't include baboon_id label
  strip.background = element_blank(), #facet laabel background white, not grey
  panel.spacing = unit(0, "lines") #no space between the facets aka they stack closer on top of each other
)

Fig4A_panel1 <- ggplot(m2, aes(x=age, y=pc1_bc, group = baboon_id2, fill = baboon_id2)) +
  geom_area(fill = "#33a02c") +
  xlab("Age (years)") +
  ylab("Individual Bray-Curtis PC1") +
  facet_wrap(~baboon_id2, ncol = 1, strip.position="left", scales = "free_y") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+ #adds x axis line on each panel
  xlim(min(m2$age)-0.1,max(m2$age)+0.05) + #set xlims to just above and below min/ max age
  scale_y_continuous(breaks = c(0)) + 
  scale_x_continuous(expand=c(0,0), breaks = c(5,10,15,20,25)) +
  t3

Fig4A_panel2 <- ggplot(m2, aes(x=age, y=Christensenellaceae), group = baboon_id2, fill = baboon_id2) +
  geom_area(fill = "#6a3d9a") +
  xlab("Age (years)") +
  ylab("Individual Christensenellaceae CLR") +
  facet_wrap(~baboon_id2, ncol = 1, strip.position="left", scales = "free_y") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+ #adds x axis line on each panel
  xlim(1.1,22) + #set xlims to just above aand below min/ max
  scale_y_continuous(breaks = c(0)) + #quants
  scale_x_continuous(expand=c(0,0), breaks = c(5,10,15,20,25)) +
  t3

Fig4A <- plot_grid(Fig4A_panel1, Fig4A_panel2, align = "hv", ncol = 2)

#### Fig. 4B ####

#upload baboon subsample data
b1 <- read.csv("git_h2_model_nsamples_per_baboon_output.csv", row.names = 1)
b1$n_samples_per_subject <- factor(b1$n_samples_per_subject)

#calculate the mean h2 and mean percent of taxa that are heritable per subsample depth
b2 <- b1 %>%
  group_by(n_samples_per_subject,subsample_rep, h2_status, .drop=FALSE) %>%
  dplyr::summarise(percent_h2 = n(), mean_subsample_rep_h2 = mean(heritability)) %>%
  filter(h2_status == "heritable") %>%
  group_by(n_samples_per_subject) %>%
  summarise(h2_taxa_percent = mean(percent_h2), mean_h2 = mean(mean_subsample_rep_h2, na.rm = TRUE)) %>%
  mutate(host_species = "baboon") %>%
  dplyr::select(host_species, n_samples_per_subject, mean_h2, h2_taxa_percent)
b2 <- data.frame(b2)
b2 <- b2 %>%
  dplyr::select(host_species, n_samples_per_subject, mean_h2, h2_taxa_percent)

#upload human data
h1 <- read.csv("git_human_heritability_studies_all_traits.csv", row.names = 1)

#calculate the mean h2 and mean percent of taxa that are heritable per human study
h2 <- h1 %>%
  group_by(dataset, is_heritable) %>%
  dplyr::summarise(h2_traits = n(), mean_h2 = mean(h2_estimate))  %>%
  dplyr::mutate(h2_taxa_percent = 100*(h2_traits / sum(h2_traits))) %>%
  filter(is_heritable == "heritable") %>%
  mutate(host_species = "human", n_samples_per_subject = 1) 
h2 <- data.frame(h2)
h2 <- h2 %>%
  dplyr::select(host_species, n_samples_per_subject, mean_h2, h2_taxa_percent)

m1 <- rbind(h2, b2)
class(m1$n_samples_per_subject)
m1$n_samples_per_subject <- as.numeric(as.character(m1$n_samples_per_subject))
m1 <- m1 %>%
  mutate(`Host species` = host_species, `Mean heritability` = mean_h2)

t3<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line.x = element_line(size=.4), 
  axis.line.y = element_line(size=.4), 
  legend.position= c(0.4,0.3), #left, height
  legend.title=element_text(size=18), 
  legend.text=element_text(size=16),
  axis.title = element_text(color="black", size=20),
  axis.text = element_text(color="black", size=18)
)

Fig4B <- ggplot(m1, aes(x=n_samples_per_subject, y=h2_taxa_percent)) +
  geom_point(aes(fill=`Host species`, size = `Mean heritability`), colour = "black", shape = 21, alpha = 0.5) +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  scale_size_continuous(range = c(4,10), breaks = c(0.1,0.3,0.5)) +
  xlab("Samples per individual") +
  ylab("Heritable taxa (%)") +
  scale_x_continuous(breaks = c(1,2,5,10,20)) +
  t3


#### Fig. 4C ####

mods <- read.csv("git_h2_model_sampling_depth_output.csv")
mods <- mods %>%
  dplyr::filter(phenotype %in% c("asv_shannon_h","pc1_bc","g_Christensenellaceae_R-7_group","g_Prevotella_2")) %>%
  dplyr::group_by(phenotype, n_samples) %>%
  dplyr::summarise(quant100 = quantile(heritability, 1),
                   quant0 = quantile(heritability, 0),
                   quant95 = quantile(heritability, 0.95),
                   quant5 = quantile(heritability, 0.05),
                   quant75 = quantile(heritability, 0.75),
                   quant25 = quantile(heritability, 0.25)
  )

#make n samples numeric from a factor so they evenly space
mods$n_samples2 <- as.factor(mods$n_samples)
mods$n_samples2 <- as.numeric(mods$n_samples2)

t1<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line.x = element_line(size=.4), 
  axis.line.y = element_line(size=.4), 
  legend.position="top",
  legend.text = element_text(size = 16),
  axis.title = element_text(color="black", size=20),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
  axis.text.y = element_text(size = 16),
  strip.text.x = element_text(size = 16)
)

#manually label x axis
xlabs <- unique(mods$n_samples)
xloc <- unique(mods$n_samples2)

tax.labs = c("ASV Shannon's H", "Christensenellaceae R-7 grp", "Prevotella 2","Bray-Curtis PC1")

names(tax.labs) <- sort(unique(as.character(mods$phenotype)))

#add columns for levels
mods$outer <- "full data"
mods$middle <- "90th percentile"
mods$inner <- "50th percentile"


Fig4C <- ggplot(mods, aes(n_samples2)) + 
  #outer ribbon
  geom_ribbon(aes(ymin = quant0, ymax = quant100, fill = outer)) +
  #inner ribbon
  geom_ribbon(aes(ymin = quant5, ymax = quant95, fill = middle)) +
  #super inner ribbon?
  geom_ribbon(aes(ymin = quant25, ymax = quant75, fill = inner)) +
  scale_fill_manual(values=c("#005a32", "#33a02c", "#b2df8a"), name="") +
  xlab("Number of samples") +
  labs(y = expression(paste("Heritability (", h^2, ")"))) +
  facet_wrap(~phenotype,ncol=2, labeller = labeller(phenotype = tax.labs)) +
  scale_x_continuous(labels = xlabs, breaks = xloc) +
  t1 


#### Fig. 4D ####

mods <- read.csv("git_h2_model_sampling_depth_output.csv", row.names =1)
#percent models improved by pedigree. dplyr drops 0 count rows
sigs_bc <- mods %>%
  group_by(phenotype, n_samples, heritability_status) %>%
  dplyr::summarise(percent_sig = n()) %>%
  filter(heritability_status == "Heritable")

#make a column for heritable in the full dataset
sigs_full <- mods %>%
  filter(n_samples == 16234) %>%
  filter(subsample_rep == 1) %>%
  dplyr::select(phenotype, full_model_heritability_status = heritability_status)

mod2 <- merge(sigs_bc, sigs_full, by="phenotype")

#add a column for taxon or phenotype
drops <- c("asv_richness", "asv_shannon_h", "pc1_bc", "pc2_bc", "pc3_bc", "pc4_bc", "pc5_bc")
mod2$Response_variable <- ifelse(mod2$phenotype %in% drops, "Community phenotype", "Single-taxon phenotype")

#manually label x axis
#make n samples numeric from a factor so they evenly space
mod2$n_samples2 <- as.factor(mod2$n_samples)

t3<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  # panel.grid.major.y = element_line(color = "grey", size = 0.4), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line.x = element_line(size=.4), 
  axis.line.y = element_line(size=.4), 
  legend.position="top",
  legend.title=element_text(size=16), 
  legend.text=element_text(size=14),
  axis.title = element_text(color="black", size=20),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
  axis.text.y = element_text(color="black", size=16),
  strip.background = element_blank(), #facet laabel background white, not grey
)

palA <-c("#d7191c", "#2b83ba")

Fig4D <- ggplot(mod2, aes(x=n_samples2, y=percent_sig)) + 
  geom_line(aes(colour = factor(full_model_heritability_status), group = phenotype)) +
  scale_color_manual(values = palA, guide = guide_legend(title = "Full dataset",ncol=2)) +
  xlab("Number of samples") +
  ylab("Significantly heritable traits (%)") +
  facet_wrap(~Response_variable, ncol = 1, strip.position="left") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+ #adds x axis line on each panel
  t3