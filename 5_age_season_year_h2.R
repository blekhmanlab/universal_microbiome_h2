############################################################################
#### Script 5: Year, season, and host age modify heritability estimates ####
############################################################################


############################################################################
#### Section 5.1: Heritability per year ####
############################################################################

library(broom)
library(dplyr)
library(lmerTest)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggridges)
library(ggrepel) 

#### Fig 3A ####

m1 <- read.csv("git_h2_model_135yearly_output.csv", row.names = 1)
#add if h2 is significant
m1 <- m1 %>%
  mutate(h2_significance = case_when(p_adjust < 0.1 ~ "Heritable", 
                                     p_adjust >= 0.1 ~ "Not heritable"))

#add on full model h2 column
#read in heritable taxa, limit to the 15 most heritable collapsed phenotypes
mods <- read.csv("git_h2_model_283single_taxon_and_7community_phenotype_output.csv", row.names = 1)
mods <- mods %>%
  arrange(phenotype)
mods <- subset(mods, phenotype %in% unique(m1$phenotype))
mods <- mods %>%
  arrange(-heritability) %>%
  dplyr::select(phenotype, `Full dataset h2` = heritability, `Full dataset h2 SE` = heritability_se)

m1 <- merge(m1, mods, by = "phenotype")
#sort by hydro_year
m1 <- m1 %>%
  arrange(desc(-as.numeric(hydro_year)))
m1$hydro_year <- as.factor(m1$hydro_year)

#Completely clear all lines except axis lines and make background white + add a legend
t1<-theme(plot.margin = unit(c(0, 0, 0.1, 0), "cm"),                               
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(size=.4), 
          axis.line.y = element_line(size=.4), 
          legend.position="top",
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(color="black", size=20),
          strip.background = element_blank(), #facet laabel background white, not grey
          strip.text.x = element_text(size = 12)
)

tax.labs = c("Prevotellaceae",
             "CAG-873",
             "Candidatus Methanogranum",
             "Christensenellaceae R-7 group",
             "Collinsella",
             "Family XIII AD3011 group",
             "Libanicoccus",
             "Prevotella 2",
             "Prevotella 9",
             "Rikenellaceae RC9 gut group",
             "Ruminococcaceae UCG-009",
             "Ruminococcaceae UCG-011",
             "Mollicutes RF39",
             "WCHB1-41",
             "Bray-Curtis PC1"
)


names(tax.labs) <- sort(unique(as.character(m1$phenotype)))

#colors
palA <-c("#d7191c", "#2b83ba")
palB <- c("#A9A9A9")

#add a dummy column for grey SE bars
m1$dummy_grey <- "Heritability +/- SE across all years"


Fig3A <- ggplot(m1, aes(x=as.numeric(hydro_year))) + 
  geom_ribbon(aes(ymin=`Full dataset h2`-`Full dataset h2 SE`, ymax=`Full dataset h2` + `Full dataset h2 SE`, fill = dummy_grey)) + 
  scale_fill_manual(values = palB, guide = guide_legend(title = "")) + 
  geom_line(aes(y=heritability, group = phenotype), size = 0.25, color = "black") +
  geom_point(aes(y=heritability, colour=h2_significance), size = 2.5) +
  scale_colour_manual(values = palA, guide = guide_legend(title = "Heritability per year", ncol = 1)) +
  geom_errorbar(aes(y=heritability, ymin=heritability-heritability_se, ymax=heritability+heritability_se, colour=h2_significance), width=.4, position=position_dodge(.9)) + 
  scale_x_continuous(breaks = 1:9, labels = levels(m1$hydro_year)) +
  facet_wrap(~phenotype, ncol=3, labeller = labeller(phenotype = tax.labs), scales = 'free_y') +
  xlab("Hydrological year") +
  labs(y = expression(paste("Heritability (", h^2, ") +/- SE"))) +
  t1


############################################################################
#### Section 5.2: Heritability per season ####
############################################################################

#### Season analyses ####

## Test if wet and dry season h2 estimates are correlated for all phenotypes 
m1 <- read.csv("git_h2_model_100wetseason_output.csv", row.names = 1)
m2 <- read.csv("git_h2_model_100dryseason_output.csv", row.names = 1)
m2 <- merge(m1, m2, by = "phenotype")
cor.test(m2$heritability.x, m2$heritability.y)

## Limit to phenotypes heritable in both seasons
m3 <- m2 %>%
  filter(p_adjust.x < 0.1) %>%
  filter(p_adjust.y < 0.1)
cor.test(m3$heritability.x, m3$heritability.y)

## Test if wet season h2 < dry season h2
t.test(m3$heritability.x, m3$heritability.y, paired = TRUE)

## Test if wet season Vp > dry season Vp for all phenotypes
m2 <- m2 %>%
  mutate(Vp.x = VA.x+VI.x+Vmom.x+Vplate.x+VR.x) %>%
  mutate(Vp.y = VA.y+VI.y+Vmom.y+Vplate.y+VR.y) %>%
  #exclude richness and shannon's h, as their Vps are several orders of magnitude higher
  filter(phenotype != "asv_shannon_h") %>%
  filter(phenotype != "asv_richness")
t.test(m2$Vp.x, m2$Vp.y, paired = TRUE)


#### Fig 3B ####

m1 <- read.csv("git_h2_model_100wetseason_output.csv", row.names = 1)
m1 <- m1 %>%
  mutate(h2_wet = heritability, h2_se_wet = heritability_se, p_adjust_wet = p_adjust)

m2 <- read.csv("git_h2_model_100dryseason_output.csv", row.names = 1)
m2 <- m2 %>%
  mutate(h2_dry = heritability, h2_se_dry = heritability_se, p_adjust_dry = p_adjust)

mods <- merge(m1, m2, by = "phenotype")

#make 1 pedigree status column
mods$Heritable <- ifelse((mods$p_adjust_dry < 0.1) & (mods$p_adjust_wet < 0.1), "Both seasons", 
                         ifelse((mods$p_adjust_dry < 0.1) & (mods$p_adjust_wet >= 0.1), "Dry season only",
                                ifelse((mods$p_adjust_dry >= 0.1) & (mods$p_adjust_wet < 0.1), "Wet season only",
                                       ifelse((mods$p_adjust_dry >= 0.1) & (mods$p_adjust_wet >= 0.1), "Neither season", "error"))))

#order levels
mods$Heritable <- as.factor(mods$Heritable)
mods$Heritable <- factor(mods$Heritable, levels(mods$Heritable)[c(1,2,4,3)])


#Completely clear all lines except axis lines and make background white + add a legend
t3<-theme(plot.margin = unit(c(0, 0, 0.1, 0), "cm"),                               
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(size=.4), 
          axis.line.y = element_line(size=.4), 
          legend.position=c(0.025, 0.9),
          legend.text = element_text(size = 14),
          axis.title = element_text(color="black", size=20),
          axis.text = element_text(size = 18)
)


#pick colors
pal<-c("#d7191c", "#fdbf6f", "#b2df8a", "#3288bd")

Fig3B <- ggplot(mods, aes(y=h2_dry, x=h2_wet)) +
  geom_abline(intercept = 0, slope = 1, colour = "darkgrey", weight = 2, linetype = 'dashed') +
  geom_errorbarh(aes(xmax = h2_wet+h2_se_wet, xmin = h2_wet-h2_se_wet, colour=Heritable), size=0.6) +
  geom_errorbar(aes(ymax = h2_dry+h2_se_dry, ymin = h2_dry-h2_se_dry, colour=Heritable), size=0.6) +
  geom_point(aes(fill=Heritable), size = 3, shape = 21) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  geom_smooth(method="lm", se=TRUE, color = "black") +
  labs(y = expression(paste("Dry season heritability (", h^2, ") +/- SE")), x = expression(paste("Wet season heritability (", h^2, ") +/- SE"))) +
  t3

#### Fig 3C ####

#limit it to only microbes heritable in both seasons
mods2 <- mods %>%
  filter(Heritable == "Both seasons")

#Note direction
mods2$Direction <- ifelse(mods2$h2_dry > mods2$h2_wet, "Dry > wet", "Wet > dry")

wet_mods <- mods2 %>%
  dplyr::select(phenotype, h2 = h2_wet, Direction) %>%
  mutate(Season = "Wet")
dry_mods <- mods2 %>%
  dplyr::select(phenotype, h2 = h2_dry, Direction) %>%
  mutate(Season = "Dry")
mods2 <- rbind(wet_mods, dry_mods)

t3<-theme(plot.margin = unit(c(0, 0, 0.1, 0), "cm"),                               
          plot.background = element_blank(), 
          panel.grid.major.x =  element_blank(), 
          panel.grid.major.y =  element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(size=.4), 
          legend.position=c(0.3, 0.9),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=20),
          axis.text = element_text(size = 18)
)

#pick the colors
linepal <- c("#fdbf6f","#b2df8a")
pointpal<-c("#d7191c", "#3288bd")

Fig3C <- ggplot(mods2, aes(y=h2, x=Season, group = phenotype)) +
  geom_line(aes(colour = Direction), size = 1) +
  geom_point(aes(fill = Direction), shape=21, size = 2) +
  labs(y = expression(paste("Heritability (", h^2, ")"))) +
  scale_colour_manual(values = linepal) +
  scale_fill_manual(values = linepal) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3), labels = c(0.0, 0.1, 0.2, 0.3), limits = c(0.0,0.3)) +
  t3


#### Fig 3D ####

m1 <- read.csv("git_metadata_and_community_phenotypes.csv")

palA <-c("#fdbf6f", "#b2df8a")

t1<-theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),                               
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(size=.4), 
          legend.position="none",
          #  axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(color="black", size=20),
          axis.text = element_text(size = 18)
)

Fig3D <- ggplot(m1, aes(y=diet_shannon_h, x=season)) +
  geom_boxplot(aes(fill = season)) + 
  scale_fill_manual(values = palA) +
  ylab("Dietary Shannon's H") +
  t1



############################################################################
#### Section 3.3: Heritability per age class ####
############################################################################

#### Testing which phenotypes change with age ####

ups <- read.csv("git_h2_model_1400age_output.csv", row.names = 1)

#make model number a factor, then numeric
ups$model_number <- as.factor(ups$model_number)
ups$model_number <- factor(ups$model_number, levels(ups$model_number)[c(1:2,7:14,3:6)])
ups$model_number_numeric <- as.numeric(ups$model_number)

#loop through and run a linear model for each taxon
dfHour = ups %>% 
  group_by(phenotype) %>%
  do(fitHour = lm(heritability ~ model_number_numeric, data = .))

dfHourCoef = tidy(dfHour, fitHour)

#limit to phenotypes with increasing h2 with age
ups2 <- dfHourCoef %>%
  filter(term == "model_number_numeric") %>%
  filter(p.value < 0.05) %>%
  filter(estimate > 0) %>%
  arrange(-estimate)

#### Testing changes in Vp, Vr, and Va with age ####

m1 <- ups %>%
  #exclude phenotypes where variance is several orders of magnitude higher
  filter(phenotype != "asv_shannon_h") %>%
  filter(phenotype != "asv_richness") %>%
  mutate(VP = VA+VI+Vmom+Vplate+VR)

lm1 <- lmer(VA ~ model_number_numeric + n_samples + n_snames + (1|phenotype), data=m1)
lm1 <- lmer(VR ~ model_number_numeric + n_samples + n_snames + (1|phenotype), data=m1)
lm1 <- lmer(VP ~ model_number_numeric + n_samples + n_snames + (1|phenotype), data=m1)

#### Testing if diet diversity changes with age ####

mods <- read.csv("git_diet_by_age.csv", row.names = 1)
mods$model_number <- as.factor(mods$age_window)
mods$model_number <- factor(mods$model_number, levels(unique(mods$model_number))[c(1:2,7:14,3:6)])

lm1 <- lmer(diet_shannon_h ~ as.numeric(model_number) + (1|baboon_id), data = subset(mods, season == "dry"))
lm1 <- lmer(diet_shannon_h ~ as.numeric(model_number) + (1|baboon_id), data = subset(mods, season == "wet"))


#### Testing if grooming partner diversity change with age ####

f1 <- read.csv("git_grooming_by_age.csv", row.names = 1)
lm1 <- lmer(grooming_partner_diversity ~ age_in_years + (1|social_groups) + (1|hydrological_years) + (1|baboon_id), data = f1)
lm1 <- lmer(grooming_partner_richness ~ age_in_years + (1|social_groups) + (1|hydrological_years) + (1|baboon_id), data = f1)

#### Testing if microbiome diversity change with age ####

m1 <- read.csv("git_metadata_and_community_phenotypes.csv", header = TRUE, row.names = 1)
p1 <- read.csv("git_pedigree.csv", header = TRUE, row.names = 1)
m1 <- merge(m1, p1, by.x = "baboon_id", by.y = "ID")

m1$collection_date <- as.Date(m1$collection_date, "%Y-%m-%d")
m1$plate <- as.factor(m1$plate)
m1$collection_date_numeric <- as.numeric(m1$collection_date)
m1$month <- as.factor(m1$month)
m1$hydro_year <- as.factor(m1$hydro_year)

lm1 <- lmer(asv_shannon_h ~ age + collection_date_numeric + readcount + sex + rain_month_mm + social_group + group_size + diet_PC1 + diet_PC2 + diet_PC3 + diet_PC4 + diet_PC5 + diet_PC6 + diet_PC7 + diet_PC8 + diet_PC9 + diet_PC10 + diet_PC11 + diet_PC12 + diet_PC13
            + post_pcr_dna_ng
            + hydro_year
            + month
            + (1|baboon_id) + (1|MOTHER) + (1|plate), data = m1)


#### Fig 3E ####

mods_ups <- merge(ups2, ups, by = "phenotype")

#reverse levels
mods_ups$model_number_rev <- mods_ups$model_number
mods_ups$model_number_rev <- factor(mods_ups$model_number_rev, levels(mods_ups$model_number_rev)[c(14,13,12,11,10,9,8,7,6,5,4,3,2,1)])
levels(mods_ups$model_number_rev)

#Completely clear all lines except axis lines and make background white
t3<-theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),                               
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(size=.4), 
          axis.line.y = element_line(size=.4), 
          legend.position="none",
          axis.title = element_text(color="black", size=20),
          axis.text = element_text(color="black", size=18)
)

#set colors
palA <- c("#a1d99b", "#31a354")

Fig3E <- ggplot(mods_ups, aes(
  x=heritability,
  y=model_number_rev,
  fill=model_number_rev,
  rel_min_height=0.000000000001
)) +
  # Adding "stat = 'density' means the bandwidth (i.e. bin size)
  # is calculated separately for each category, rather than for the
  # dataset as a whole
  stat_density_ridges(
    scale = 2,
    quantile_lines = TRUE,
    quantiles = 2
  ) +
  scale_fill_cyclical(values=palA) +
  labs(x=expression(paste("Heritability (", h^2, ")")), y='Age class (years)') +
  #add a vertival line at the median
  geom_vline(
    xintercept=median(mods_ups$heritability),
    col="yellow", linetype="dashed", size=1
  ) +
  t3

#### Fig 3F ####

#limit to phenotypes with 10 biggest changes
ups10 <- data.frame(ups2)
ups10 <- ups10 %>%
  arrange(-estimate) %>%
  slice(1:10)

mods10 <- merge(ups10, ups, by = "phenotype")
#add a column for h2
mods10$Heritable <- ifelse(mods10$p_adjust < 0.1, "Heritable", "Not heritable")

t3<-theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),                              
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_blank(), 
          axis.line.y = element_line(size=.4), 
          legend.position=c(0.025, 0.9),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          axis.title = element_text(color="black", size=20),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
          axis.text.y = element_text(size = 18),
          legend.text.align = 0 #left align the legend text
)


tax.labs1 <- mods10 %>%
  filter(model_number == "13-27")  %>%
  dplyr::select(phenotype, heritability) %>%
  distinct() %>%
  arrange(-heritability)
tax.labs1 <- droplevels(tax.labs1)

tax.labs1$phenotype

# [1] o_Mollicutes_RF39               g_Ruminococcaceae_UCG-002       g_Prevotella_2                 
# [4] g_Bifidobacterium               g_Ruminococcaceae_UCG-014       asv_shannon_h                  
# [7] g_Ruminococcaceae_NK4A214_group g_Ruminococcaceae_UCG-011       g_Senegalimassilia             
# [10] g_Helicobacter 

tax.labs = c(expression("Mollicutes RF39"),
             expression(italic("Ruminococcaceae UCG-002")),
             expression(italic("Prevotella 2")),
             expression(italic("Bifidobacterium")),
             expression(italic("Ruminococcaceae UCG-014")),
             expression("ASV Shannon's H"),
             expression(italic("Ruminococcaceae NK4A214 group")),
             expression(italic("Ruminococcaceae UCG-011")),
             expression(italic("Senegalimassilia")),
             expression(italic("Helicobacter")))

#use the brewer palette 
colourCount = length(unique(tax.labs1$phenotype))
getPalette = colorRampPalette(brewer.pal(length(unique(tax.labs1$phenotype)), "Set1"))

palA <-c("#d7191c", "#2b83ba")


Fig3F <- ggplot(mods10, aes(y=heritability, x=model_number)) +
  geom_line(aes(group=phenotype, colour = phenotype), size = 1, show.legend = FALSE) +
  scale_colour_manual(values = getPalette(colourCount)) +
  geom_text_repel(data = mods10 %>% filter(model_number == "13-27") %>% arrange(heritability), label = rev(tax.labs), aes(color = phenotype),direction = "y", hjust = 0, show.legend = FALSE, xlim = c(14.5,NA), segment.size = NA, size = 4) + 
  geom_point(aes(fill=Heritable), shape = 21, size = 2.5) +
  scale_fill_manual(values = palA, guide = guide_legend(title = "Heritability per age class", ncol = 1)) +
  labs(y = expression(paste("Heritability (", h^2, ")"))) +
  xlab("Age class (years)                                                ") + 
  coord_cartesian(xlim = c(1, 20)) + #make the x axis extend long enough that all the text fits
  geom_segment(aes(x=0, xend=14, y=-Inf, yend=-Inf)) + #but only show line til number
  t3