############################################################################
#### Script 3: Genetic effects on the gut microbiome are near-universal ####
############################################################################

library(dplyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(RColorBrewer)
library(scales)

#### Calculating per-sample prevalence vs heritability ####

t1 <- readRDS("git_PA_table.RDS")
t1$sample <- NULL
t1 <- t1 %>%
  summarise_each(funs(sum))
t1 <- as.data.frame(t(t1))
t1$percent_prevalence <- t1$V1/16234

t2 <- read.csv("git_h2_model_744presenceabsence_output.csv")
t1 <- merge(t2,t1, by.x = "phenotype", by.y = "row.names")
cor.test(t1$heritability, t1$percent_prevalence)

#### Comparing heritability estimates between transformations ####

# Single-taxon phenotype (relative abundance) vs CLR transformation
t1 <- read.csv("git_h2_model_283clr_output.csv")
t1 <- t1 %>%
  mutate(clr_h2 = heritability)
t2 <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv")
t2 <- t2 %>%
  mutate(ab_h2 = heritability)
t2 <- merge(t1, t2, by = "phenotype")

cor.test(t2$ab_h2, t2$clr_h2)
t.test(t2$ab_h2, t2$clr_h2, paired=TRUE)

# Single-taxon phenotype (relative abundance) vs presence/absence 
t1 <- read.csv("git_h2_model_744presenceabsence_output.csv")
t1 <- t1 %>%
  mutate(pa_h2 = heritability)
t2 <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv")
t2 <- t2 %>%
  mutate(ab_h2 = heritability)
t2 <- merge(t1, t2, by = "phenotype")

cor.test(t2$ab_h2, t2$pa_h2)
t.test(t2$ab_h2, t2$pa_h2, paired=TRUE)

#### Comparing additive genetic variance vs maternal and individual variances ####
t2 <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv")
#calculate maternal and individual percent variance
t2 <- t2 %>%
  mutate(mom_prop = Vmom/(VA+VI+Vplate+VR+Vmom)) %>%
  mutate(id_prop = VI/(VA+VI+Vplate+VR+Vmom))

#heritability vs maternal variance
t.test(t2$heritability, t2$mom_prop, paired = TRUE)
#heritability vs individual variance
t.test(t2$heritability, t2$id_prop, paired = TRUE)


#### Figure 2A ####

mods <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv", row.names = 1)
drops <- c("asv_richness", "asv_shannon_h", "pc1_bc", "pc2_bc", "pc3_bc", "pc4_bc", "pc5_bc")
mods$Response_variable <- ifelse(mods$phenotype %in% drops, "Community phenotype", "40 most heritable taxa")
mods_top <- mods %>%
  filter(!phenotype %in% drops) %>%
  arrange(-heritability) 

## add on ASV taxonomy info
phyla1 <- readRDS("git_ASV_table.RDS")
phyla1 <- phyla1 %>%
  dplyr::select(domain, phylum, class, order, family, genus, asv_id)
long_lab <- merge(mods_top, phyla1, by.x = "phenotype", by.y = "asv_id")

long_lab <- long_lab %>%
  mutate(domain = paste("d_", domain, sep = ""), phylum = paste("p_", phylum, sep = ""), class = paste("c_", class, sep = ""), order= paste("o_", order, sep = ""), family= paste("f_", family, sep = ""), genus= paste("g_", genus, sep = "")) %>%
  mutate(long_id = paste(domain, phylum, class, order, family, genus, sep = "; ")) %>%
  dplyr::select(phenotype, long_id)

#clean up the labels
long_lab$long_id <- gsub("[a-z]_NA", "NA", long_lab$long_id)
long_lab$long_id <- gsub("[a-z]_NA; ", "", long_lab$long_id)
long_lab$long_id <- gsub("; [a-z]_NA", "", long_lab$long_id)
long_lab$long_id <- gsub("NA; ", "", long_lab$long_id)
long_lab$long_id <- gsub("; NA", "", long_lab$long_id)
long_lab$long_id <- gsub("d_Bacteria; ", "", long_lab$long_id)
long_lab$long_id <- gsub("^.*; ", "", long_lab$long_id)
long_lab$levs <- gsub("_.*", "", long_lab$long_id)
long_lab$long_id <- gsub("[a-z]_", "", long_lab$long_id)

long_lab <- long_lab %>%
  mutate(new_name = paste(long_id, " (", phenotype, "; ", levs, ")", sep = ""))
long_lab$new_name <- gsub("asv", "ASV", long_lab$new_name)
long_lab <- long_lab %>%
  select(phenotype, new_name)

#add on label information
mods_top<- merge(mods_top, long_lab, by = "phenotype", all.x = TRUE)

# fix non-asv labels
mods_top<- mods_top %>%
  #rename g_ to taxa (g) and also taxa (asv)
  mutate(new_name = case_when(is.na(new_name) ~ as.character(paste(substring(phenotype,3,), " (",substr(phenotype,1,1),")",sep="")),
                              !is.na(new_name) ~ new_name))

mods_top$new_name <- gsub("_", " ", mods_top$new_name)


# use the top 40
mods_top <- mods_top %>%
  arrange(-heritability) %>%
  slice(1:40)

#order by model ordeer
mods_top$phenotype <- factor(mods_top$phenotype, levels = mods_top$phenotype[order(-mods_top$heritability)])

#remove extra columns
mods_top$levs <- NULL
mods_top$long_id <- NULL

#ditto community phenotypes
mods_ph <- mods %>%
  arrange(-heritability) %>%
  filter(phenotype %in% drops) %>%
  mutate(new_name = case_when(phenotype == "asv_richness" ~ "ASV richness",
                              phenotype == "asv_shannon_h" ~ "ASV Shannon's H",
                              phenotype == "pc1_bc" ~ "Bray-Curtis PC1",
                              phenotype == "pc2_bc" ~ "Bray-Curtis PC2",
                              phenotype == "pc3_bc" ~ "Bray-Curtis PC3",
                              phenotype == "pc4_bc" ~ "Bray-Curtis PC4",
                              phenotype == "pc5_bc" ~ "Bray-Curtis PC5"
  ))

#order by modsel ordeer
mods_ph$phenotype <- factor(mods_ph$new_name, levels = mods_ph$new_name[order(-mods_ph$heritability)])
levels(mods_ph$phenotype)

mods_b <- rbind(mods_top, mods_ph)
mods_b <- droplevels(mods_b)

#pick colors
palBcol <-c("#6a3d9a", "#33a02c")

t3<-theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),                          
          plot.background = element_blank(), 
          panel.grid.major.x = element_line(color = "grey", size = 0.4), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(size=.4), 
          axis.line.y = element_line(size=.4), 
          legend.position="top",
          # legend.position=c(-1,1),
          legend.justification='left',
          legend.title=element_blank(), 
          legend.text=element_text(size=14),
          axis.title.y = element_blank(),
          axis.title.x = element_text(color="black", size=20),
          axis.text = element_text(color="black", size=16)
)


Fig2A <- ggplot(mods_b, aes(x=phenotype, y=heritability)) +
  geom_point(aes(colour=Response_variable), size = 3) +
  scale_x_discrete(limits = rev(levels(mods_b$phenotype)), labels = rev(mods_b$new_name)) +
  geom_errorbar(aes(ymax = heritability+heritability_se, ymin = heritability-heritability_se, colour=Response_variable), size=0.6) +
  scale_colour_manual(values = palBcol, guide = guide_legend(title = "",ncol=1)) +
  labs(y = expression(paste("Heritability (", h^2, ") +/- SE"))) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), label = c(0, "", 0.1, "", 0.2, "")) +
  coord_flip() + 
  t3


#### Figure 2B ####

# 283 single-taxa phenotypes
rh <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv", row.names = 1)
nope <- c("asv_richness", "pc1_bc", "pc2_bc", "pc3_bc", "pc4_bc", "pc5_bc", "asv_shannon_h")
rh <- rh %>%
  filter(!phenotype %in% nope) %>%
  mutate(model_type = "Relative abundance") %>%
  mutate(model_type_significance = case_when(p_adjust < 0.1 ~ "sig", 
                                             p_adjust >= 0.1 ~ "NS")) %>%
  dplyr::select(model_type, model_type_significance, heritability)

# clr 
clr <- read.csv("git_h2_model_283clr_output.csv", row.names = 1)
clr <- clr %>%
  mutate(model_type = "CLR") %>%
  mutate(model_type_significance = case_when(p_adjust < 0.1 ~ "sig", 
                                             p_adjust >= 0.1 ~ "NS")) %>%
  dplyr::select(model_type, model_type_significance, heritability)

# philr
philr <- read.csv("git_h2_model_138philr_output.csv", row.names = 1)
philr <- philr %>%
  mutate(model_type = "PhILR") %>%
  mutate(model_type_significance = case_when(p_adjust < 0.1 ~ "sig", 
                                             p_adjust >= 0.1 ~ "NS")) %>%
  dplyr::select(model_type, model_type_significance, heritability)

# binomal
bin <- read.csv("git_h2_model_744presenceabsence_output.csv", row.names = 1)
bin <- bin %>%
  mutate(model_type = "Presence/absence") %>%
  mutate(model_type_significance = case_when((VA_constraint == "Boundary") ~ "NS", 
                                             heritability - heritability_se > 0 ~ "sig", 
                                             heritability - heritability_se < 0  ~ "NS")) %>%
  dplyr::select(model_type, model_type_significance, heritability)

#combine them
m1 <- rbind(rh, clr, philr, bin)
m1$model_type_significance2 <- paste(m1$model_type, m1$model_type_significance, sep = "_")

#color by h2 but all of them are the same
palD2 <- c('#6a3d9a', '#cab2d6')

#make sure NS listed first so it shows up in front on the graph
m1$model_type_significance <- as.factor(m1$model_type_significance)
m1$model_type_significance <- factor(m1$model_type_significance,levels(m1$model_type_significance)[c(2,1)])

#make model type a factor so it's
m1$model_type <- as.factor(m1$model_type)
m1$model_type <- factor(m1$model_type,levels(m1$model_type)[c(4,1,2,3)])

#make summary text about # heritable at p<0.1
sig_counts <- m1 %>%
  group_by(model_type, model_type_significance) %>%
  dplyr::summarise(n())
l2 <- c("273/283 (96%) heritable", "280/283 (99%) heritable", "132/138 (96%) heritable", "704/744 (95%) heritable")
labs1 <- data.frame(model_type = unique(m1$model_type))
labs1$label_it <- l2

#make a label mean
dMean2 <- m1 %>%
  filter(model_type_significance == "sig") %>%
  group_by(model_type) %>%
  summarise(MN = signif(mean(heritability), digits = 2))

t2<-theme(plot.margin = margin(6, 6, 0, 6),                              
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(size=.4), 
          axis.line.y = element_blank(), 
          legend.position="none",
          axis.title.x = element_text(color="black", size=20),
          axis.text.x = element_text(color="black", size=18),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          # strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(angle = 0, size=18), #1 is right aligned, anfle 0 is rotates
          panel.spacing = unit(0, "lines") #no space between the facets aka they stack closer on top of each other
)

Fig2B <- ggplot(m1, aes(
  x=heritability,
  fill= model_type_significance,
  colour = model_type_significance
)) + 
  geom_histogram(position="identity", binwidth = 0.007) +
  scale_color_manual(values = palD2) +
  scale_fill_manual(values = alpha(palD2, 0.8)) +
  facet_wrap( ~ model_type, ncol=1, strip.position = "top") +
  labs(x = expression(paste("Heritability (", h^2, ")"))) +
  geom_text(data = labs1, aes(x = 0.18, y = 29, label = label_it), size = 5, inherit.aes = FALSE) +
  geom_vline(data = dMean2, aes(xintercept = as.numeric(MN)), col="yellow", size=1) +
  geom_text(data = dMean2, aes(x = as.numeric(MN), y = 45, label = MN), size = 5, inherit.aes = FALSE) +
  t2


#### Figure 2C ####

# 283 single-taxa phenotypes + 7 community phenotypes
mods <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv", row.names = 1)
mods3 <- mods %>%
  mutate(Technical = 100*(Vplate/(Vplate+Vmom+VA+VI+VR)),
         Maternal = 100*(Vmom/(Vplate+Vmom+VA+VI+VR)),
         `Additive genetic` = 100*(VA/(Vplate+Vmom+VA+VI+VR)),
         `Individual identity` = 100*(VI/(Vplate+Vmom+VA+VI+VR))
  ) %>%
  dplyr::select(phenotype, Technical, Maternal, `Additive genetic`, `Individual identity`)

#melt the data to get the Components as a column
mods3 <- reshape2::melt(mods3)
colnames(mods3) <- c("phenotype", "Component", "est")

#add on the p value column
mods2 <- subset(mods, select = c("phenotype", "p_adjust"))
mods3 <- merge(mods2, mods3, by.x = "phenotype")
colnames(mods3) <- c("Model", "pedigree_p", "component", "est")

plotdat <- mods3

# Need to reorder levels so they plot in the right order left to right: technical, maternal, baboon id, additive genetic
plotdat$component <- as.factor(plotdat$component)
plotdat$component <- factor(plotdat$component,levels(plotdat$component)[c(1,2,4,3)])

#add a column for if it's taa or not
not_tax <- c("pc1_bc", "pc2_bc", "pc3_bc", "pc4_bc", "pc5_bc", "asv_richness", "asv_shannon_h")

#reorder phenotypes so it's by highest heritability, and group by taxa vs community phenotype
h_order <- mods
h_order$status <- ifelse(h_order$phenotype %in% not_tax, "Non-taxa phenotype", "Taxa")
h_order$tax_level <- substr(h_order$phenotype, 0,2)
#add in level, order levels
h_order <- h_order %>%
  mutate(tax_level2 = case_when(status == "Non-taxa phenotype" ~ "composition",
                                tax_level == "p_" ~ "phylum",
                                tax_level == "c_" ~ "class",
                                tax_level == "o_" ~ "order",
                                tax_level == "f_" ~ "family",
                                tax_level == "g_" ~ "genus",
                                tax_level == "as" ~ "ASV")) %>%
  mutate(tax_level_numeric = case_when(status == "Non-taxa phenotype" ~ 7,
                                       tax_level == "p_" ~ 6,
                                       tax_level == "c_" ~ 5,
                                       tax_level == "o_" ~ 4,
                                       tax_level == "f_" ~ 3,
                                       tax_level == "g_" ~ 2,
                                       tax_level == "as" ~ 1)) %>%
  arrange(tax_level_numeric, heritability) %>%
  group_by(status) %>%
  dplyr::mutate(variable_order = 1:n()) %>%
  dplyr::select(phenotype, variable_order, status, tax_level2, heritability) 

plotdat <- merge(plotdat, h_order, by.x = "Model", by.y = "phenotype")


pal <- c('#a1dab4','#016450','#41b6c4','#225ea8')

## no legend
th2 <- theme(
  axis.line.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  legend.position = "none",
  axis.text.y = element_blank(),
  axis.title.x = element_text(color="black", size=20),
  axis.text.x = element_text(color="black", size=18),
  axis.ticks.y = element_blank(),
  strip.background = element_blank() #facet laabel background white, not grey
  , strip.text.y =  element_blank() #remove teh facet labels
  , strip.text.x = element_text()
)

#community phenotype plot
pC_1 <- ggplot(subset(plotdat, status == "Non-taxa phenotype"), aes(x = variable_order, y = est, fill=component)) +
  geom_bar(stat = "identity",colour="white", size=0.25) +
  scale_fill_manual(values = pal) +
  guides(fill = guide_legend(reverse = TRUE, title = "Variance component")) +  #make additive genetic top
  coord_flip(clip = "off") + #wont' clip text labels
  labs(y="Percent variance explained", x=" ") +
  annotate("text", x = 10, y = 30,label = "Community phenotype (n=7)" , color="black", size=5) +
  ##add in dummy lines so everything lines up
  annotate("segment", x = 2, xend = 4, y = -1, yend = -1, colour = "white") + 
  annotate("text", x = 3, y = -3, size = 3, angle = 90, label = "b", color = "white") + 
  th2


line_info <- h_order %>%
  group_by(tax_level2) %>%
  summarise(line_min = min(variable_order), line_max = max(variable_order), line_mid = mean(variable_order))

#taxa plot + legend
th3 <- theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
             #axis.line.x = element_blank(), 
             axis.line.y = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank(),
             legend.title = element_text(size=14),
             legend.text = element_text(size=12),
             legend.position = c(0.47, 0.075),
             axis.text.y = element_blank(),
             axis.title.y = element_text(color="black", size=20),
             axis.title.x = element_text(color="black", size=20),
             axis.text.x = element_text(color="black", size=18),
             axis.ticks.y = element_blank(),
             strip.background = element_blank() #facet laabel background white, not grey
             , strip.text.y =  element_blank() #remove teh facet labels
             , strip.text.x = element_text()
)

pC_2 <- ggplot(subset(plotdat, status != "Non-taxa phenotype"), aes(x = variable_order, y = est, fill=component)) +
  geom_bar(stat = "identity",colour="white", size=0.25) +
  scale_fill_manual(values = pal) +
  guides(fill = guide_legend(reverse = TRUE, title = "Variance component")) +  #make additive genetic top
  coord_flip(clip = "off") + #wont' clip text labels
  labs(y="Percent variance explained", x="Microbial phenotype") +
  annotate("text", x = 290, y = 15,label = "Single-taxon phenotype (n=283)" , color="black", size=5) +
  scale_x_discrete(expand=c(0,0)) +
  
  ## phylum ## 272-283
  annotate("segment", x = 272, xend = 283, y = -1, yend = -1, colour = "darkgrey") + 
  annotate("text", x = 280, y = -3, size = 5, angle = 90, label = "Phylum") +
  
  ## class ## 253 - 271
  annotate("segment", x = 253, xend = 271, y = -1, yend = -1, colour = "darkgrey") + 
  annotate("text", x = 258, y = -3, size = 5, angle = 90, label = "Class") +
  
  ## order ## 232 - 252
  annotate("segment", x = 232, xend = 252, y = -1, yend = -1, colour = "darkgrey") + 
  annotate("text", x = 238, y = -3, size = 5, angle = 90, label = "Order") +
  
  ## family ## 201 - 231
  annotate("segment", x = 201, xend = 231, y = -1, yend = -1, colour = "darkgrey") + 
  annotate("text", x = 216, y = -3, size = 5, angle = 90, label = "Family") +
  
  ## genus ## 140 - 200
  annotate("segment", x = 140, xend = 200, y = -1, yend = -1, colour = "darkgrey") + 
  annotate("text", x = 170, y = -3, size = 5, angle = 90, label = "Genus") +
  
  ## asv ## 253 - 271
  annotate("segment", x = 1, xend = 139, y = -1, yend = -1, colour = "darkgrey") + 
  annotate("text", x = 70, y = -3, size = 5, angle = 90, label = "ASV") +
  
  th3

Fig2C <- plot_grid(pC_1, pC_2, align = "hv", ncol = 1, rel_heights = c(1.5,10))