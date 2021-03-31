
############################################################################
#### Script 2: Running the heritability models ####
############################################################################

# NOTE: ASReml-R is available to purchase from VSN International. This section describes how to run the animal model if you have accees to ASReml-R. If you don't, skip to section 3 for interpreting model outputs. 

library(pedantics)
library(asreml)
library(dplyr)
library(nadiv)

#### 2.1 Running the animal model on 7 community phenotypes ####

#upload metadata
m1 <- read.csv("git_metadata_and_community_phenotypes.csv", header = TRUE, row.names = 1)

#asreml requires dates as numeric
m1$collection_date <- as.Date(m1$collection_date, "%Y-%m-%d")
m1$collection_date_numeric <- as.numeric(m1$collection_date)

#make plate, month, and hydrologic year factors
m1$hydro_year <- as.factor(m1$hydro_year)
m1$month <- as.factor(m1$month)
m1$plate <- as.factor(m1$plate)

#read in the pedigree file, one column for "id" (child), one for "dam" (mother), one for "sire" (father); NA for unknown values
ped = read.csv("git_pedigree.csv", header = TRUE, row.names = 1)

#prepare the pedigree
ped2 <- fixPedigree(ped)
#note that ped2 outputs column names as "id, dam, sire" but the dam and sire labels are switched for the MOTHER and FATHER columns. asreml requires column order to be ID, FATHER, MOTHER
colnames(ped2) <-  c("ID", "FATHER", "MOTHER")

#create the inverse of the relationship matrix from the pedigree file
ainv<-asreml.Ainverse(ped2)$ginv

#add the mother and father columns to the metadata
m1 <- merge(m1, ped2, by.x = "baboon_id", by.y = "ID") 


#upload the microbiome phenotype file (in this case, compositional phenotypes are in the metadata file)
tax <- m1 %>%
  dplyr::select(sample, asv_richness, asv_shannon_h, pc1_bc, pc2_bc, pc3_bc, pc4_bc, pc5_bc)

#make columns for relevant model outputs; the variance components, heritability estimates, and heritability likelikhood ratio test chi value and p-value
t1 <- as.data.frame(matrix(nrow = ncol(tax), ncol = 14))
rownames(t1) <- colnames(tax)
colnames(t1) <- c("Vplate", "Vplate_constraint", "Vmom", "Vmom_constraint", "VA", "VA_constraint", "VI", "VI_constraint", "VR", "VR_constraint", "heritability", "heritability_se", "VA_LRT_chi", "VA_LRT_p")

# run a separate model on every taxon

for (i in 2:length(tax)){
  #make a new dataframe as sample name and 1 taxon
  a <- tax[,c(1,i)]
  #add taxon to the metadata
  a2 <- merge(m1, a, by = "sample")
  #call furthest right column (the taxon) 'taxa'
  a2$taxa <- a2[,c(39)]
  
  modely<-asreml(fixed=taxa~ 1 + collection_date_numeric + readcount + age + sex + rain_month_mm + social_group + group_size + diet_PC1 + diet_PC2 + diet_PC3 + diet_PC4 + diet_PC5 + diet_PC6 + diet_PC7 + diet_PC8 + diet_PC9 + diet_PC10 + diet_PC11 + diet_PC12 + diet_PC13
                 + post_pcr_dna_ng
                 + hydro_year
                 + month
                 , random=  ~ped(baboon_id,var=T,init=1)+ plate + MOTHER + ide(baboon_id,var=T,init=1)
                 , data=a2
                 , ginverse=list(baboon_id=ainv)
                 , na.method.X="include", na.method.Y="omit")
  
  
  #example model output from summary(modely)$varcomp
  #                                    gamma    component    std.error   z.ratio constraint
  # plate!plate.var             1.359187e+00 4.337450e+03 4.395693e+02  9.867501   Positive
  # MOTHER!MOTHER.var           1.615056e-07 5.153980e-04 5.867053e-06 87.846142   Boundary
  # ped(baboon_id, var = T)!ped 2.279129e-01 7.273177e+02 1.290650e+02  5.635281   Positive
  # ide(baboon_id, var = T)!id  1.128601e-01 3.601601e+02 7.975380e+01  4.515899   Positive
  # R!variance                  1.000000e+00 3.191209e+03 3.632725e+01 87.846142   Positive
  
  ### constraint = 'Boundary' means the variance component contributes ~0 to the model
  ### constraint = 'Positive' means the variance component makes some contribution to the model
  
  #calculate heritability and heritability SE from variance components as VA/(Vplate + Vmom + VA + VI + VR)
  modely_h2 <- nadiv:::pin(modely, h2 ~ V3 / ( V1 + V2 + V3 + V4 + V5))
  
  #add the relevant model output to a table
  t1[i,1] <- summary(modely)$varcomp[c("plate!plate.var"),c("component")]
  t1[i,2] <- as.character(summary(modely)$varcomp[c("plate!plate.var"),c("constraint")])
  t1[i,3] <- summary(modely)$varcomp[c("MOTHER!MOTHER.var"),c("component")]
  t1[i,4] <- as.character(summary(modely)$varcomp[c("MOTHER!MOTHER.var"),c("constraint")])
  t1[i,5] <- summary(modely)$varcomp[c("ped(baboon_id, var = T)!ped"),c("component")]
  t1[i,6] <- as.character(summary(modely)$varcomp[c("ped(baboon_id, var = T)!ped"),c("constraint")])
  t1[i,7] <- summary(modely)$varcomp[c("ide(baboon_id, var = T)!id"),c("component")]
  t1[i,8] <- as.character(summary(modely)$varcomp[c("ide(baboon_id, var = T)!id"),c("constraint")])
  t1[i,9] <- summary(modely)$varcomp[c("R!variance"),c("component")]
  t1[i,10] <- as.character(summary(modely)$varcomp[c("R!variance"),c("constraint")])
  t1[i,11] <- modely_h2[1,1] #h2
  t1[i,12] <- modely_h2[1,2] #h2_se
  
  #run the model without pedigree
  modely_no_ped<-asreml(fixed=taxa~ 1 + collection_date_numeric + readcount + age + sex + rain_month_mm + social_group + group_size + diet_PC1 + diet_PC2 + diet_PC3 + diet_PC4 + diet_PC5 + diet_PC6 + diet_PC7 + diet_PC8 + diet_PC9 + diet_PC10 + diet_PC11 + diet_PC12 + diet_PC13
                        + post_pcr_dna_ng
                        + hydro_year
                        + month
                        , random= ~ plate + MOTHER + ide(baboon_id,var=T,init=1)
                        , data=a2
                        , ginverse=list(baboon_id=ainv)
                        , na.method.X="include", na.method.Y="omit")
  
  #compare models using a likelihood ratio test
  t1[i,13] <- 2*(modely$loglik-modely_no_ped$loglik)  #VA lrt
  t1[i,14] <- 1-pchisq(2*(modely$loglik-modely_no_ped$loglik),1) #VA p
}

#remove the top, blank row
t_community = t1[-1,]

#### 2.2 Running the animal model on 283 single-taxon phenotypes (relative abundance) found in >50% of samples ####

# Same methods as (2.1), except the 1. response variable is modified to fixed=asin(sqrt(taxa)) (see below) and 2. the tax table is the relative abundance table

tax <- readRDS("git_relative_abundance_table.RDS")

modely<-asreml(fixed=asin(sqrt(taxa))~ 1 + collection_date_numeric + readcount + age + sex + rain_month_mm + social_group + group_size + diet_PC1 + diet_PC2 + diet_PC3 + diet_PC4 + diet_PC5 + diet_PC6 + diet_PC7 + diet_PC8 + diet_PC9 + diet_PC10 + diet_PC11 + diet_PC12 + diet_PC13
               + post_pcr_dna_ng
               + hydro_year
               + month
               , random=  ~ped(baboon_id,var=T,init=1)+ plate + MOTHER + ide(baboon_id,var=T,init=1)
               , data=a2
               , ginverse=list(baboon_id=ainv)
               , na.method.X="include", na.method.Y="omit")


#rbind the community and single-taxon phenotypes output
t1 <- rbind(t_community, t_single_taxon)

# Correct p-values for multiple testing using a Benjamini-Hochberg correction at a 10% FDR

#For all taxa with Positive VA constraints, run a 10% FDR correction using Benjamini-Hochberg (BH)
high_h <- subset(t1, (VA_constraint != "Boundary"))
high_h$p_adjust <- p.adjust(high_h$VA_LRT_p, method = "BH")
#add on the adjusted p
high_h <- subset(high_h, select = c("p_adjust"))
t1 <- merge(t1, high_h, by = "row.names", all.x = TRUE)
#assign taxa with Boundary VA constraint a p-value of 1
t1$p_adjust[is.na(t1$p_adjust)] <- 1

#save a copy of the output (this file is provided)
write.csv(t1, "git_h2_model_283single_taxon_and_7community_phenotype_output.csv")


#### 2.3 Running the animal model on CLR-transformed 283 single-taxon phenotypes found in >50% of samples ####

# Same methods as (2.1), except the taxa table is the CLR table
tax <- readRDS("git_CLR_table.RDS")

#save a copy of the output (this file is provided)
write.csv(t1, "git_h2_model_283clr_output.csv")

#### 2.4 Running the animal model on PhILR-transformed 138 ASVs found in >50% of samples ####

# Same methods as (2.1), except the tax table is the philr table
tax <- readRDS("git_philr_table.RDS")

#save a copy of the output (this file is provided)
write.csv(t1, "git_h2_model_138philr_output.csv")


#### 2.5 Running the animal model on 744 presence/absence phenotypes ####

# Note that: 1. model family and maxiter are added, 2. heritability is calculated differently, and (3) likelihood ratio tests cannot be used to determine model significance so an alternate approach is used

# Prep the metadata and pedigreefollowing (2.1)

tax <- readRDS("git_PA_table.RDS")

# make columns for relevant model outputs
t1 <- as.data.frame(matrix(nrow = ncol(tax), ncol = 10))
rownames(t1) <- colnames(tax)
colnames(t1) <- c("Vplate", "Vplate_constraint", "Vmom", "Vmom_constraint", "VA", "VA_constraint", "VI", "VI_constraint", "heritability", "heritability_se")

# run a separate model on every taxon in the input file

for (i in 2:length(tax)){
  #make a new dataframe as sample name and 1 taxon
  a <- tax[,c(1,i)]
  #add taxon to the metadata
  a2 <- merge(m1, a, by = "sample")
  #call furthest right column (the taxon) 'taxa'
  a2$taxa <- a2[,c(39)]
  
  modely<-asreml(fixed=taxa~ 1 + collection_date_numeric + readcount + age + sex + rain_month_mm + social_group + group_size + diet_PC1 + diet_PC2 + diet_PC3 + diet_PC4 + diet_PC5 + diet_PC6 + diet_PC7 + diet_PC8 + diet_PC9 + diet_PC10 + diet_PC11 + diet_PC12 + diet_PC13
                 + post_pcr_dna_ng
                 + hydro_year
                 + month
                 , random=  ~ped(baboon_id,var=T,init=1)+ plate + MOTHER + ide(baboon_id,var=T,init=1)
                 , family = asreml.binomial(link = "logit")
                 , data=a2
                 , ginverse=list(baboon_id=ainv)
                 , na.method.X="include", na.method.Y="omit"
                 , maxiter = 20)
  
  #example model output from summary(modely)$varcomp
  
  #                                gamma    component  std.error  z.ratio constraint
  # MOTHER!MOTHER.var       9.139588e-08 9.139588e-08         NA       NA   Boundary
  # plate!plate.var         1.914451e-01 1.914451e-01 0.02948733 6.492454   Positive
  # ped(sname, var = T)!ped 2.387441e-01 2.387441e-01 0.05831224 4.094236   Positive
  # ide(sname, var = T)!id  6.965575e-02 6.965575e-02 0.04067438 1.712522   Positive
  # R!variance              1.000000e+00 1.000000e+00         NA       NA      Fixed
  
  #calculate heritability and heritability SE from variance components. Note that residual variance (R!variance) is a dummy output (constraint = Fixed), and the variance associated with the link function (logit link variance = 3.29) must be used in the denominator instead
  
  modely_h2 <- nadiv:::pin(modely, h2 ~ V3 / ( V1 + V2 + V3 + V4 + 3.29))
  
  t1[i,1] <- summary(modely)$varcomp[c("plate!plate.var"),c("component")]
  t1[i,2] <- as.character(summary(modely)$varcomp[c("plate!plate.var"),c("constraint")])
  t1[i,3] <- summary(modely)$varcomp[c("MOTHER!MOTHER.var"),c("component")]
  t1[i,4] <- as.character(summary(modely)$varcomp[c("MOTHER!MOTHER.var"),c("constraint")])
  t1[i,5] <- summary(modely)$varcomp[c("ped(baboon_id, var = T)!ped"),c("component")]
  t1[i,6] <- as.character(summary(modely)$varcomp[c("ped(baboon_id, var = T)!ped"),c("constraint")])
  t1[i,7] <- summary(modely)$varcomp[c("ide(baboon_id, var = T)!id"),c("component")]
  t1[i,8] <- as.character(summary(modely)$varcomp[c("ide(baboon_id, var = T)!id"),c("constraint")])
  t1[i,9] <- modely_h2[1,1] #h2
  t1[i,10]<- modely_h2[1,2] #h2_se
}

#remove the top, blank row
t1 = t1[-1,]

#calculate significance as if additive genetic variance (VA) makes a non-zero contribution to the model, and if the heritability standard error bars do not cross 0 
t1 <- t1 %>%
  mutate(model_type_significance = case_when((VA_constraint == "Boundary") ~ "NS", 
                                             heritability - heritability_se > 0 ~ "sig", 
                                             heritability - heritability_se < 0  ~ "NS")) 

#save a copy of the output (this file is provided)
write.csv(t1, "git_h2_model_744presenceabsence_output.csv")