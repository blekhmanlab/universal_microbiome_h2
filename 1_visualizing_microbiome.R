
############################################################################
#### Script 1: Visualizing the microbiome data set  ####
############################################################################

library(dplyr)
library(ggplot2)
library(RColorBrewer)

#### Figure 1A ####

#order by who has the youngest sample
m1 <- read.csv("git_metadata_and_community_phenotypes.csv", header = TRUE, row.names = 1)

# order by youngest age
m_count <- m1 %>%
  group_by(baboon_id) %>%
  summarise(n_samples = min(age)) %>%
  arrange(n_samples) %>%
  ungroup()
m_count$sample_count_order <- 1:nrow(m_count)

m2 <- merge(m1, m_count, by = "baboon_id")
m2$sex <- ifelse(m2$sex == "M", "male", "female")

#set the theme and plot
t3<-theme(plot.margin = margin(0, 6, 0, 6),                               
          plot.background = element_blank(), 
          panel.grid.major.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(size=.4), 
          axis.line.y = element_line(size=.4),
          axis.ticks.length.y = unit(0, "pt"),
          legend.position = c(0.6,0.9),
          legend.title = element_blank(),
          legend.text = element_text(color="black", size=20),
          axis.text.x = element_text(color="black", size=18),
          axis.title.x = element_text(color="black", size=20),
          axis.title.y = element_text(color="black", size=20)
)

Fig1A <- ggplot(m2, aes(x=age, y=-sample_count_order)) +
  geom_point(size = 1, shape=16, alpha = 0.5, aes(colour = sex)) +
  scale_colour_manual(values = c("#4DAF4A", "#000075")) + #green navy
  guides(colour = guide_legend(override.aes = list(size=3))) +
  xlab("Age (years)") +
  ylab("Baboon individual") +
  scale_x_continuous(expand=c(0,0), breaks = c(5,10,15,20,25)) +
  scale_y_discrete(expand=c(0.01,0)) +
  t3

#### Figure 1D ####

m1 <- read.csv("git_nonhumanprimate_and_human_microbiome_phyla.csv", row.names = 1)
#note that the name Kiritimatiellaeota replaces Verrucomicrobia in Mann and Berer studies, as these are the same phylum

# make taxa a factor
m1$taxa <- as.factor(m1$taxa)

## group by Phyl Group
phyl_order <- m1 %>%
  group_by(monkey) %>%
  sample_n(1) %>%
  mutate(Phy2 = case_when(PhylGroup == "Lemur" ~ 1,
                          PhylGroup == "New_World" ~ 2,
                          PhylGroup == "Ape" ~ 3,
                          PhylGroup == "Old_World" ~ 4)) %>%
  arrange(Phy2, monkey) %>%
  dplyr::select(Phy2, PhylGroup, monkey) 

phyl_order$monkey_order <- 1:nrow(phyl_order)

m1 <- merge(m1, phyl_order, by = c("PhylGroup", "monkey"))

#order by primate group
m1$monkey <- as.factor(m1$monkey)
m1$monkey2 <- factor(m1$monkey, levels = unique(m1$monkey[order(m1$monkey_order)]))
m1$CommonName <- as.factor(m1$CommonName)
m1$CommonName2 <- factor(m1$CommonName, levels = unique(m1$CommonName[order(m1$monkey_order)]))

#use the same color palette as Fig 1E
phy_cols <- c('#3cb44b', '#ffe119', '#4363d8', '#f58231', 
              '#911eb4', '#42d4f4', '#f032e6', '#bfef45', 
              '#fabebe', '#469990')


#theme and plot
t2<-theme(                             
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.y = element_line(size=.4), 
  legend.position="none",
  legend.title = element_blank(),
  axis.title.y = element_text(color="black", size=16),
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 14),
  axis.text.y = element_text(size = 14),
  axis.ticks.length.x = unit(0, "pt")
)

Fig1D <- ggplot(m1, aes(x=CommonName2, y=mean_phylum_abundance, fill=taxa)) + 
  geom_bar(stat="identity") +
  ylab("Mean relative abundance \n of phylum in microbiome") +
  xlab("") +
  scale_fill_manual(values = phy_cols) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) +
  scale_x_discrete(position = "bottom") +
  ##lemurs ##
  annotate("segment", x = 1, xend = 2, y = 104, yend = 104) + 
  annotate("text", x = 1.5, y = 109, label = "Lemur") +
  ## new world ##
  annotate("segment", x = 3, xend = 7, y = 104, yend = 104) + 
  annotate("text", x = 4.5, y = 109, label = "Platyrrhini") +
  ## apes ##
  annotate("segment", x = 8, xend = 10, y = 104, yend = 104) + 
  annotate("text", x = 8.5, y = 109, label = "Hominidae") +
  ## old world ##
  annotate("segment", x = 11, xend = 16, y = 104, yend = 104) + 
  annotate("text", x = 13.5, y = 109, label = "Cercopithecidae") +
  t2

#### Figure 1E ####

m1 <- read.csv("git_metadata_and_community_phenotypes.csv", header = TRUE, row.names = 1)

#upload ASV table
phyla1 <- readRDS("git_ASV_table.RDS")
#switch out the . for -
colnames(phyla1) <- gsub("\\.", "-", colnames(phyla1))
#limit it to phyla
phyla <- as.data.frame(phyla1[,1:7])
phyla <- subset(phyla, select = c("phylum"))
#turn phyla into a character instead of factor
phyla$phylum <- as.character(phyla$phylum)

#some taxa were unassigned. replace  blanks w that
phyla[is.na(phyla)] <- "Unassigned"
#also remove white spaces; replace with _
phyla$phylum <-gsub(" ", "_", phyla$phylum, fixed = TRUE)
#ditto parentheses
phyla$phylum <-gsub("(", "_", phyla$phylum, fixed = TRUE)
phyla$phylum <-gsub(")", "_", phyla$phylum, fixed = TRUE)

#subset to samples of interest
phyla_b <- phyla1[colnames(phyla1) %in% m1$sample]
#add taxonomy back on
phyla <- cbind(phyla,phyla_b)

#aggregate by phylum
phyla2 <- aggregate(.~phylum, data = phyla, FUN=sum)

#rescale to relative abundance, thence to percentage out of 100
phyla2[-1] <- 100*sweep(phyla2[-1],2,colSums(phyla2[-1]),`/`)

#remove unassigned row
phyla2 <- droplevels(phyla2)

#make a list of the most relatively abundant phyla
phyla_sum <- phyla2
phyla_sum$mean_abundance <- rowMeans(phyla_sum[-1])
phyla_sum2 <- subset(phyla_sum, select = c("phylum","mean_abundance"))
phyla_sum2 <- phyla_sum2[with(phyla_sum2, order(-mean_abundance)),]
#make a table to plot
phyla_sum_graph <- phyla_sum2
#make phylum the rownames
rownames(phyla_sum_graph) <- phyla_sum_graph$phylum
phyla_sum_graph$phylum <- NULL

#Keep the phyla that comprise > 0.5% of reads on average
smol <- subset(phyla_sum2, mean_abundance < 0.5)
#replace all taxa in rare list w unassigned
phyla2$abundant_phyla <- ifelse(phyla2$phylum %in% smol$phylum, "Rare or unassigned", phyla2$phylum)
#make unassisnged also into rare or unassigned
phyla2$abundant_phyla <- ifelse(phyla2$abundant_phyla == "Unassigned", "Rare or unassigned", phyla2$abundant_phyla)
#remove the other phylum column
phyla2$phylum <- NULL
#aggregate again so rare and unassigned is 1 row
phyla3 <- aggregate(.~abundant_phyla, data = phyla2, FUN=sum)

#melt it
phyla22 <- reshape2::melt(phyla3, id.vars = "abundant_phyla")

#tack on a columns for categorical variables of interest
phyla23 <- merge(phyla22, m1, by.x= "variable", by.y = "sample")

#order by date
phyla24 <- subset(phyla23, select = c("variable", "collection_date"),abundant_phyla == "Firmicutes")
phyla24 <- phyla24[with(phyla24, order(collection_date)),]
phyla24$variable_order <- 1:nrow(phyla24)
phyla24 <- subset(phyla24, select = c("variable", "variable_order"))
phyla23 <- merge(phyla23, phyla24, by = "variable")

phyla23 <- phyla23[with(phyla23, order(variable_order)),]
phyla23$legend_order <- 1:nrow(phyla23)

#add dashed bars for years
yrs <- phyla23 %>%
  dplyr::select(year = hydro_year, variable_order) %>%
  distinct()

#X location of year label
vars <- yrs %>% 
  group_by(year) %>%
  summarise(median_var = (median(c(max(variable_order),min(variable_order)))))
yrs_x <- vars$median_var

#Y location of year label
yrs_y <- rep(-7,14)
#yr labels
yrs_labs <- unique(yrs$year)

#create year label text
yrs_labs2 <- c("", 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012,2013)

#location for years bar
yrs_bar <- yrs %>%
  group_by(year) %>%
  arrange(year) %>%
  slice(1L)

#graph it

#use color palette from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
phy_cols <- c('#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabebe', '#469990')

#set the theme + plot
t1<-theme(plot.margin = margin(0, 0, 6, 6),                               
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size=.4), 
          legend.position="right",
          legend.title = element_blank(),
          legend.text=element_text(size=18),
          axis.title.x = element_text(color="black", size=20),
          axis.title.y = element_text(color="black", size=20),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color="black", size=18)
)

Fig1E <- ggplot(phyla23, aes(x=factor(variable_order), y=value, fill=abundant_phyla)) + 
  geom_bar(stat="identity") +
  ylab("Relative abundance of \n phylum in microbiome") +
  xlab("Collection date") +
  t1 +
  scale_fill_manual(values = phy_cols) +
  geom_vline(xintercept = yrs_bar$variable_order, linetype = 2) +
  annotate("text", x = yrs_x+200, y = yrs_y, label = yrs_labs2, size=5, angle = 90, vjust = 0) +
  coord_cartesian(clip="off")