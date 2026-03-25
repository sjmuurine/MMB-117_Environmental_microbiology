
####################################################################################

#   .................................  MMB-117 .....................................#
#
###################################################################################

# libraries and packages

library(reshape2)
library(dplyr)
library(tidyr)
library(MASS)
library(vegan)
library(ggplot2)
library(ggthemes)
library(colorspace)
library(colorRamps)
library(RColorBrewer)
library(gridExtra)
library(fitdistrplus)
library(logspline)
library(car)
library(lme4) 
#library(multcomp)
library(scales)
library(viridis)
library(ggrepel)
library(pheatmap)
#library(psych)

###################################################################################

#                             Housekeeping                                        #

###################################################################################


#R-environment is cleared

rm(list=ls())


setwd("/Users/jmuurine/Desktop/Students")

# read in matadata


MMB117metadata <- read.table("MMB117_Sample_sheet_2025.txt", header=TRUE, stringsAsFactors = TRUE, sep="\t")

View(MMB117metadata)
MMB117metadata<-MMB117metadata[1:25,1:20]

# We can live with that although there are some funny things

# read in data

ASV_table <- readRDS("seqtab.rds")

str(ASV_table)
View(ASV_table)

tax_table <- readRDS("taxa.rds")
str(tax_table)

View(tax_table)

identical(rownames(tax_table), colnames(ASV_table))
#[1] TRUE

head(rownames(tax_table))

head(colnames(ASV_table))

# Let's change the column and row names

colnames(ASV_table) <- paste0("ASV", seq(ncol(ASV_table)))

rownames(tax_table) <- paste0("ASV", seq(nrow(tax_table)))

# I will store the row sums for later use if I need:
ASV_counts_per_sampleDF <- as.data.frame(rowSums(ASV_table))

# I will check if there are ASVs that are detected only less than 3 times

dim(ASV_table)
# [1]    26 20262

dim(ASV_table[, which(colSums(ASV_table)<3)])
#[1]  26 606

# One could remove these:
#ASV_table_over3<-ASV_table[, which(colSums(ASV_table)>3)]

str(ASV_table)
class(ASV_table)


# TSS normalization to get realtive abundances for data exploration

ASV_tableRA.mat<- apply(ASV_table, 1, function(i) i/sum(i))
ASV_tableRA.mat<-t(ASV_tableRA.mat)

View(ASV_tableRA.mat)

# I will remove the negative control from sequencing lab

ASV_tableRA.mat<-ASV_tableRA.mat[2:nrow(ASV_tableRA.mat),]
View(MMB117metadata)

# the Sample names are different...
# I will change them.

MMB117metadata$Sample
#[1] GA_01        GA_02        GA_03        GA_04        GA_05        GA_06        PA_01        PA_02        PA_03        PA_04        PA_05        PA_06       
#[13] FI_01        FI_02        FI_03        FI_04        FI_05        FI_06        FO_01        FO_02        FO_03        FO_04        FO_05        FO_06       
#[25] Neg. control

rownames(ASV_tableRA.mat)
#[1] "CTRL-TATCCAGT" "FI-01"         "FI-02"         "FI-03"         "FI-04"         "FI-05"         "FI-06"         "FO-01"         "FO-02"         "FO-03"        
#[11] "FO-04"         "FO-05"         "FO-06"         "GA-01"         "GA-02"         "GA-03"         "GA-04"         "GA-05"         "GA-06"         "PA-01"        
#[21] "PA-02"         "PA-03"         "PA-04"         "PA-05"         "PA-06"        

MMB117Samples<- c("Neg. control", "FI_01", "FI_02", "FI_03",  "FI_04",  "FI_05",  "FI_06",  
                  "FO_01" , "FO_02", "FO_03", "FO_04",  "FO_05" , "FO_06",
                  "GA_01",  "GA_02" ,  "GA_03",  "GA_04",  "GA_05", "GA_06",
                  "PA_01",  "PA_02" , "PA_03" ,"PA_04"  ,"PA_05",   "PA_06")
length(MMB117Samples)
#[1] 25

rownames(ASV_tableRA.mat)<-MMB117Samples

# a basic NMDS plot usually tells if there is something funny with the data

plot(metaMDS(ASV_tableRA.mat), type="text", display="sites")

#Warning message:
#In metaMDS(ASV_tableRA.mat) :
 # stress is (nearly) zero: you may have insufficient data

# because of the negative control, as the plot shows.

plot(metaMDS(ASV_tableRA.mat[2:25,]), type="text", display="sites")

# looks pretty good! Separated by sites, so we can see that the microbiomes are different

# For the future I want to put the samples in the same order in the metadata and ASV data.

ASV_tableRA.mat_o<-ASV_tableRA.mat[match(MMB117metadata$Sample,row.names(ASV_tableRA.mat)),]

# Now I can also get Shannons diversity index and put that into metadata

MMB117metadata$ASV_divSha<-diversity(ASV_tableRA.mat_o)




###################################################################################

#                                 Data exploration                                #

###################################################################################

# According to
# Zuur, A.F., Ieno, E.N. and Elphick, C.S. (2010), 
# A protocol for data exploration to avoid common statistical problems. 
# Methods in Ecology and Evolution, 1: 3-14. https://doi.org/10.1111/j.2041-210X.2009.00001.x

#############################################################################################


# I want to start from ***INDEPENDENCE***, since it is the most important assumption:

# We had only one sampling point.
# Thus we don't really have "traditional" independence problem which would arise from repeated measurements.
# If we would have taken samples before and after some treatment from the same 
# sites, we would have to take this into account in our statistical analyses
# (use mixed models, etc.)

############################################################################################

# --- Outliers in X & Y--- #

# X = metadata and things that are usually on x-axis

# Y = the numerical data we use, in this case we have  microbiome data that are either counts (non-normalized), 
# relative abundance (TSS-normalized), or clr-transformed normally distributed numbers that can be used 
# in statistical analyses

# metadata first.
str(MMB117metadata)

#'data.frame':	25 obs. of  20 variables:
# $ Sample                  : Factor w/ 26 levels "","FI_01","FI_02",..: 14 15 16 17 18 19 21 22 23 24 ...
#$ Site                    : Factor w/ 6 levels "","-","Field",..: 5 5 5 5 5 5 6 6 6 6 ...
#$ SnowDepth               : Factor w/ 6 levels "","0 cm","10 cm",..: 3 3 1 1 2 1 4 4 1 3 ...
#$ WetWeight               : num  10 10 11.2 17.4 10 ...
#$ DryWeight               : num  6.51 6.47 6.84 13.51 5.57 ...
#$ X80C_DryWeight          : num  NA NA 6.79 13.43 5.44 ...
#$ DNAWeight               : num  306 237 501 285 340 ...
#$ Moisture                : num  58 58 65.1 32 84.7 ...
#$ pH_W                    : num  6.61 5.76 6.18 6.37 6.19 6.37 5.8 6.37 6.6 6.4 ...
#$ pH_Ca                   : num  5.63 5.24 5.72 5.77 5.58 5.62 5.34 5.57 5.81 5.74 ...
#$ SOM                     : num  13.5 12 7.6 3.8 17.6 ...
#$ ResponsiblePerson       : Factor w/ 10 levels "","Arttu","Erin",..: 4 4 6 6 5 5 9 9 9 3 ...
#$ Notes                   : Factor w/ 21 levels "","clay soil with too much water",..: 7 12 14 20 6 19 3 1 9 13 ...
#$ X10_5_TSA_Plate         : Factor w/ 17 levels "","0","1","103",..: 10 4 12 5 14 17 15 11 6 13 ...
#$ X10_6_TSA_Plate         : int  1 12 1 1 3 5 6 2 4 10 ...
#$ X10_7_TSA_Plate         : int  0 3 3 0 0 1 0 0 0 1 ...
#$ X10_4_Malt_Plate        : int  18 30 11 7 17 12 1 2 0 8 ...
#$ X10_5_Malt_Plate        : int  8 6 1 1 3 1 0 0 0 2 ...
#$ X10_6_Malt_Plate        : int  0 1 0 0 1 0 0 0 0 1 ...
#$ Other_plate_observations: Factor w/ 6 levels "","10-5 malt plate had fungal growth",..: 1 1 1 1 2 5 1 1 1 1 ...
#$ ASV_divSha              : num  7.19 7.33 7.2 7.07 7.27 ...


# So lets plot all variables that make sense to look at 

ggplot(MMB117metadata, aes(Sample, SnowDepth, label=Site)) +
 geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Snow depth")

# not recorded for all, maybe can be left put from the analysis

ggplot(MMB117metadata, aes(Sample,  WetWeight, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Wet weight (g)")

#Warning messages:
#1: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
#2: Removed 1 row containing missing values or values outside the scale range (`geom_text_repel()`). 

# negative control


ggplot(MMB117metadata, aes(Sample,  DryWeight, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Dry weight (air dried) (g)")

#Warning messages:
#1: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
#2: Removed 1 row containing missing values or values outside the scale range (`geom_text_repel()`). 

# negative control


ggplot(MMB117metadata, aes(Sample,  X80C_DryWeight, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Owen dry weight, pore water dried out (g)")

#Warning messages:
#1: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
#2: Removed 1 row containing missing values or values outside the scale range (`geom_text_repel()`). 

# negative control


ggplot(MMB117metadata, aes(Sample,  DNAWeight, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
ggtitle("What was weighed for DNA extraction (mg)")

#Warning messages:
#1: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
#2: Removed 1 row containing missing values or values outside the scale range (`geom_text_repel()`). 

# negative control


ggplot(MMB117metadata, aes(Sample,  Moisture, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Soil water % (?)")

#Warning messages:
#1: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
#2: Removed 1 row containing missing values or values outside the scale range (`geom_text_repel()`). 

# negative control

#Gravimetric water content would be better than "moisture", you could change the column. 
#There are some that are very wet. 
#Interesting to see if their microbiome is completely different because of the snow for instance


ggplot(MMB117metadata, aes(Sample,  pH_W, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Water-extraxted pH")


#Warning messages:
#1: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
#2: Removed 1 row containing missing values or values outside the scale range (`geom_text_repel()`). 

# negative control


ggplot(MMB117metadata, aes(Sample,  pH_Ca, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("CaCl2-extraxted pH")

# nice grouping!!

# With both pH's actually

#Warning messages:
#1: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 
#2: Removed 1 row containing missing values or values outside the scale range (`geom_text_repel()`). 

# negative control


ggplot(MMB117metadata, aes(Sample,  SOM, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Soil organic matter (LOI)")

# few samples that are under or above but not too extreme I would say
# some missing data as well


ggplot(MMB117metadata, aes(Sample,  ASV_divSha, label=Site)) +
  geom_point(color="black", size=3)+
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust=0, vjust=0.5))+
  geom_text_repel(color="black",size=4,max.overlaps = Inf )+
  ggtitle("Shannon's diversity index, Y data")

# negative control is an outlier, as it should

# The NMDS plot also shows if there are otliers in the samples (in terms of Y data)
plot(metaMDS(ASV_tableRA.mat[2:25,]), type="text", display="sites")
# and not really

plot(metaMDS(ASV_tableRA.mat[2:25,]))

# this shows both the species and sites and from this we cans see that different 
# species associate to different sites

              ##----------------------------- ##

# Good! looks like there are no weird extreme cases in the data.

###################################################################

# --- Homogenity of variances --- #

# With statistical analyses we usually compare means, so if the black lines are in different levels,
# that is fine. 
# but the variability of the observations should be similar. 
# i.e. the size of the boxes should be about the same

# So lets check this with metadata and Shannon's diversity index

ggplot(subset(MMB117metadata, Sample!="Neg. control"), aes(Site,  Moisture)) + 
  geom_boxplot() +
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black"))+
  theme(legend.position="none")


# whoops! Good thing we are not studying change of soil water % in response to land use!
# but I wouldn't build a linear model using soil water % as an explaining variable.
# non-parametric statistical test would be probably ok.


ggplot(subset(MMB117metadata, Sample!="Neg. control"), aes(Site,  pH_W)) + 
  geom_boxplot() +
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black"))+
  theme(legend.position="none")

# same here but less extreme

ggplot(subset(MMB117metadata, Sample!="Neg. control"), aes(Site,  pH_Ca)) + 
  geom_boxplot() +
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black"))+
  theme(legend.position="none")

# this could work with liner model as an X-variable

ggplot(subset(MMB117metadata, Sample!="Neg. control"), aes(Site,  SOM)) + 
  geom_boxplot() +
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black"))+
  theme(legend.position="none")

#Warning message:
#Removed 3 rows containing non-finite outside the scale range (`stat_boxplot()`). 

# Some samples dont have SOM -values

# Gas station is "in another planet", so SOM as an X-variable of a linear model would probably be a bad idea.

ggplot(subset(MMB117metadata, Sample!="Neg. control"), aes(Site,  ASV_divSha)) + 
  geom_boxplot() +
  theme_bw(base_size = 14)+
  theme(panel.border = element_rect(colour = "black"))+
  theme(legend.position="none")


# Also diversity shows heterogeneity in the variances, but this probably would work in linear models as an Y-variable

                          ##----------------------------- ##

# We can take these into account by selecting statistical models that that do not require homogeneity.
# If we would really need to use linear models, we could also transform these variables
# or delete observations that are outliers and go with for instance 3 observations per site.
# This is the reason why even 6 observations (6 samples) per site can bee too little!


###############################################################################################

# --- Normality --- #

#Often people say that the data should be normally distributed, 
#if you want to apply linear models.  
#This is not exactly true: The normality should be met at each X value. 
#This is of difficult to check if the data does not have multiple 
#observations for each X and we don't. 
#In other words, the Y-data contains the effects of all the explanatory variables (all metadata).
#Therefore it is misleading to assess normality according the Y data.
#Better choice is to apply a model and plot the pooled residuals. 
#Since the linear models should be LINEAR, the residuals should have 
#the same distance to fitted values (predicted means). If they do not, 
#then the model is not linear and the plot of residuals against fitted 
#values will show a curve or some other pattern. This means that 
#the RESIDUALS should be normally distributed, if linear models are used.

#We know this about normal distribution:
# "It is a continuous probability distribution in which the random variable can take any value." 
#With TSS-normalization, our values are "relative abundances", so something between 0 and 1.
#So relative abundances don't meet the normality assumption and we cannot use a model that 
#assumes normality with TSS-normalized data.

# Relative abundances are proportions, so the probability distribution to use with them would be 
# Beta distribution, but according to my experience, Gamma distribution works also 
# Gamma distribution is used to model positive continuous variables that represent time intervals between events, 
# the Beta distribution is used to model continuous variables that represent proportions or probabilities
# Gamma models work very well, I have had trouble with Beta

#######################################################################################################
# But what should be normally distributed even though it would be calculated from relative abundances:
# Shannon's diversity index
#####################################################################################################
#####################################################################################################
# Some metadata variables (that are not %) would be normally distributed, for instance pH.
#####################################################################################################

# tehäänkö tähän liittyen jotain? äänestän et ei. Tää voi mennä teorialla, menee liian vaikeeksi muuten ja ohi asian

# ---  "Fixed X" --- #

#This assumption means that you know the values of X for each sample in advance. 
#You can have e.g. different fertilization rates, treatments, etc. 
#This assumption is violated, if you for instance study some parameter of animals,
#which happen to be found somewhere and your X is for example the age of animals. In other words, 
#this assumption is violated if the X is randomly selected.

# we don't have this problem

###########################################################################################################

# ---  Zero trouble in  Y-data  --- #

# if most observations are zeros

# we probable have many zeros, at least in the ASV and genera levels

range(ASV_tableRA.mat_o)
# [1] 0.0000000 0.1498224

plot(c(ASV_tableRA.mat_o))

# and it seems we do. So if we would like to apply linear models using taxa as the response variables,
# we might need to consider zero-inflated models

###########################################################################################################

# --- Collinearity in X and # --- Relationships Y & X --- #

#If X-observations correlate, for instance weight and length, or water depth and distance to the shoreline. 
#If collinearity is ignored, one is likely to end up with a confusing statistical analysis in which nothing 
#is significant, but where dropping one covariate can make the others significant, 
#or even change the sign of estimated parameters.

# collinearity can be detected with pairwise scatter plots comparing covariates,


plot(MMB117metadata[,c(8:11)])

# pH water and pH CaCl2 seem to correlate, as they should. But because of the collinearity issue both of them cannot be used in a model

# relationships: Y-values should be plotted against X-values to see which X-values should be included in the model. 

plot(MMB117metadata[1:24,c(8:11,21)])

# This was a lazy persons solution. We can see that these is a relationship with the diversity and the X-variables. 
# And this also makes bilogically sense.

# So if we would do linear models and use diversity as the response variable, we might need to consider
# water %, pH (one of them) and SOM %

###########################################################################################################

# --- Interactions? --- #
#

# This means that if we would do linear models, some of the X-variables might have an interaction and this can be taken into account

# Let's try to visualize if we have something like this.

# I will categorize water-content

summary(MMB117metadata$Moisture)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 29.24   38.50   42.16   47.19   54.97   84.68       1 

Water_FirstQ<-38.50
Water_Median<-42.16 
Water_ThirdQ<- 54.97

MMB117metadata$MoistureCAT<- cut(MMB117metadata$Moisture,c(0, 45,  84.68  ))

MMB117metadata$MoistureCAT  

levels(MMB117metadata$MoistureCAT )  
#[1] "(0,45]"    "(45,84.7]"

levels(MMB117metadata$MoistureCAT )  <- c("0 - 45 %","45 - 84.7 %")

MMB117metadata$MoistureCAT


ggplot(subset(MMB117metadata, Sample!="Neg. control"), aes(x=pH_Ca, y=ASV_divSha, colour=Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_color_manual(values=c("springgreen4", "chocolate4","purple4","royalblue"))+
  facet_grid(Site ~ MoistureCAT) +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  theme(panel.border = element_rect(colour = "black"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  #guides(color=guide_legend(title="Sample type"))+
  theme(legend.position="none")

# Looks like in our case we don't have enough observations for this interaction

#################################################################################


# Outcomes of the data exploration:

# No reason to remove any values (outliers)
# Heteroscedasticity in all variables except pH_Ca ans Shannon diversity index this is within acceptable levels
# Normality: If we use relative abundances, they are not normally distributed. Shannon's diversity index would be
# Fixed X: we don't have this problem
# Collinearity in X and  Relationships Y & X: We need to select only one pH and looks like there are relationships 
# between Y-data and all X-values. So perhaps we should model the influence of pH_Ca, SOM and GWC on communities.
# Interactions: we don't have enough observations for a model with interactions.

# About GWC: we can include it in the analysis for practice. But actually there is an issue and it is the
# amount of snow. We didn't have a system for removing the snow or taking the amount of snow into account.
# So GWC might not be that interesting although it could explain a lot "numerically".


#################################################################################

# Biostatistical analyses:

# One of the original research questions was to see how the diversity & community structure 
# is influenced by human activities (different sites)
# hypothesis was that the diversity would be lower in the gas station

# Based on data exploration, we also should ask what is the influence of pH_Ca, SOM and GWC on communities.

# We can answer these questions with ordination analyses
# We can model pH_Ca and diversity as linear vectors 
# and SOM % and GWC % as nonlinear surfaces (GAM)

# test that NMDS works and check the dimensions
plot(metaMDS(ASV_tableRA.mat_o[1:24,]), type="text", display="sites")

#Save the NMDS
MMB117_NMDS<-metaMDS(ASV_tableRA.mat_o[1:24,])

# add colors to metadata
levels(MMB117metadata$Site)
MMB117metadata$color<-rep(1, nrow(MMB117metadata))  
MMB117metadata<- within(MMB117metadata, color[Site=="Field"]<-"greenyellow")
MMB117metadata<- within(MMB117metadata, color[Site=="Forest"]<-"forestgreen")
MMB117metadata<- within(MMB117metadata, color[Site=="Gas_station"]<-"darkslategray4")
MMB117metadata<- within(MMB117metadata, color[Site=="Park"]<-"darkkhaki")

# I want to get rid of the negative control

MMB117metadata_noneg<-subset(MMB117metadata,Sample!="Neg. control")
MMB117metadata_noneg<-droplevels(MMB117metadata_noneg)

# environmental fitting of the NMDS with diversity and pH 

MMB117EF<-envfit(MMB117_NMDS ~ ASV_divSha + pH_Ca,MMB117metadata_noneg, permutations=999)

MMB117EF

#***VECTORS

#           NMDS1    NMDS2     r2 Pr(>r)    
#ASV_divSha -0.98618 -0.16566 0.8160  0.001 ***
 # pH_Ca      -0.96879  0.24788 0.9187  0.001 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999

# This means that pH_Ca (alone) explains 91 % of the variance 
# and that diversity (alone) explains 82 % of the variance 

# including pH_Ca and div as linear vectors and surface fitting SOM %
NMDSoplot<-ordiplot(MMB117_NMDS,type="n", xlim=c(-1.8,1.8),
                  ylim=c(-1.8,1.8),cex.axis = 1.5, cex.lab = 1.5)
with(MMB117metadata_noneg, points(MMB117_NMDS$points,pch=15, cex=2, col=MMB117metadata_noneg$color))
plot(MMB117EF, col="cyan4", cex=1)
with(MMB117metadata_noneg, ordisurf(MMB117_NMDS,SOM, add = TRUE, col = "grey20", cex=2))
#with(MMB117metadata_noneg, ordisurf(MMB117_NMDS,Moisture, add = TRUE, col = "royalblue"))
identify(NMDSoplot, "sites", labels = MMB117metadata_noneg$Sample, cex=1)
with(MMB117metadata_noneg, legend(1.6,1.8, legend= levels(Site), cex=1, bty= "n", 
       col=c("greenyellow", "forestgreen", "darkslategray4","darkkhaki"), pch=c(15,15,15,15)))
legend(1.65,0.98,"SOM %",cex=1,lty=1,col="grey20",bty= "n")

# Identify samples by clicking them, then go to console and hit "esc". 
#I would click only those that are "outliers", to get the idea how
# SOM influences on communities
# Check that your cursor is in the console if nothing happens.

# statistical interpretation: 

NMDS_SOM_surf <- ordisurf(MMB117_NMDS ~ SOM, MMB117metadata_noneg)

summary(NMDS_SOM_surf) 

#Family: gaussian 
#Link function: identity 

#Formula:
# y ~ s(x1, x2, k = 10, bs = "tp", fx = FALSE)

#Parametric coefficients:
#             Estimate  Std. Error t value Pr(>|t|)    
#(Intercept)  12.0270     0.3347   35.93 2.89e-14 ***
# ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
#  edf Ref.df     F  p-value    
#s(x1,x2) 7.164      9 8.282 0.000152 ***
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.788   Deviance explained = 86.4%
#-REML = 49.569  Scale est. = 2.3531    n = 21

# ordisurf fits a GAM model and accepts nonlinear variables
# We have an intercept model and the estimate is a mean of our response variable (SOM %)
# R-sq.(adj) =  0.788 suggests that ~79 % of the variance is explained by SOM % (pretty good!)
# Deviance explained = 86.4% indicates the goodness of fit, which is also very good in our case.

# So what the model suggests is that SOM % explains ~79 % of the variance 
# in the community structure (=beta diversity) (p < 0.05 )


# Interestingly, the diversity is growing to the direction of the gas station!
# Why?
# I is lowest in the forest and field sites
# Why??

# same plot with GWC

# including pH_Ca and div as linear vectors and surface fitting GWC %
NMDSoplot<-ordiplot(MMB117_NMDS,type="n", xlim=c(-1.8,1.8),
                    ylim=c(-1.8,1.8),cex.axis = 1.5, cex.lab = 1.5)
with(MMB117metadata_noneg, points(MMB117_NMDS$points,pch=15, cex=2, col=MMB117metadata_noneg$color))
plot(MMB117EF, col="cyan4", cex=1)
#with(MMB117metadata_noneg, ordisurf(MMB117_NMDS,SOM, add = TRUE, col = "grey20", cex=2))
with(MMB117metadata_noneg, ordisurf(MMB117_NMDS,Moisture, add = TRUE, col = "royalblue", cex=2))
identify(NMDSoplot, "sites", labels = MMB117metadata_noneg$Sample, cex=1)
with(MMB117metadata_noneg, legend(1.6,1.8, legend= levels(Site), cex=1, bty= "n", 
                                  col=c("greenyellow", "forestgreen", "darkslategray4","darkkhaki"), pch=c(15,15,15,15)))
legend(1.65,0.98,"GWC %",cex=1,lty=1,col="royalblue",bty= "n")

# Identify samples by clicking them. I would click only those that are "outliers", to get the idea how 
# GWC influences on communities
# hit "esc" when you are ready. Check that your cursor is in the terminal if nothing happens.


# do the statistical interpretation for GWC yourself: 

NMDS_GWC_surf <- ordisurf(MMB117_NMDS ~ Moisture, MMB117metadata_noneg)

summary(NMDS_GWC_surf) 

#Family: gaussian 
#Link function: identity 

#Formula:
 # y ~ s(x1, x2, k = 10, bs = "tp", fx = FALSE)

#Parametric coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   47.191      1.559   30.28 2.86e-16 ***
 # ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
 # edf Ref.df    F  p-value    
#s(x1,x2) 5.941      9 5.26 0.000263 ***
 # ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.673   Deviance explained = 75.7%
#-REML = 89.604  Scale est. = 58.306    n = 24



# We can also model influence of the variables on metadata with Permanova:

adonis2(ASV_tableRA.mat_o[1:24,] ~ Site + pH_Ca + SOM + Moisture, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = ASV_tableRA.mat_o[1:24, ] ~ Site + pH_Ca + SOM + Moisture, data = MMB117metadata_noneg, permutations = 9999, by = "terms", na.action = na.omit)
#Df SumOfSqs      R2      F Pr(>F)    
#Site      3   3.7750 0.57535 8.0953 0.0001 ***
#pH_Ca     1   0.1677 0.02556 1.0789 0.3460    
#SOM       1   0.2339 0.03565 1.5049 0.1425    
#Moisture  1   0.2084 0.03176 1.3406 0.1932    
#Residual 14   2.1762 0.33167                  
#Total    20   6.5612 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Permanova suggests that only Site is significant. 
# Lets see what happens if we leave the site away.

adonis2(ASV_tableRA.mat_o[1:24,] ~ pH_Ca + SOM + Moisture, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = ASV_tableRA.mat_o[1:24, ] ~ pH_Ca + SOM + Moisture, data = MMB117metadata_noneg, permutations = 9999, by = "terms", na.action = na.omit)
#Df SumOfSqs      R2      F Pr(>F)    
#pH_Ca     1   2.0646 0.31467 9.8115 0.0001 ***
#SOM       1   0.3296 0.05024 1.5664 0.1213    
#Moisture  1   0.5896 0.08986 2.8018 0.0089 ** 
#Residual 17   3.5773 0.54523                  
#Total    20   6.5612 1.00000                  
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now SOM is not significant! Let's see what happens if we leave away GWC (moisture).


adonis2(ASV_tableRA.mat_o[1:24,] ~ pH_Ca + SOM, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = ASV_tableRA.mat_o[1:24, ] ~ pH_Ca + SOM, data = MMB117metadata_noneg, permutations = 9999, by = "terms", na.action = na.omit)
#Df SumOfSqs      R2      F Pr(>F)    
#pH_Ca     1   2.0646 0.31467 8.9187 0.0001 ***
#SOM       1   0.3296 0.05024 1.4239 0.1606    
#Residual 18   4.1669 0.63509                  
#Total    20   6.5612 1.00000                  
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# so according to permanova SOM is not significant.

# what happens if we leave out pH?

adonis2(ASV_tableRA.mat_o[1:24,] ~  SOM, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = ASV_tableRA.mat_o[1:24, ] ~ SOM, data = MMB117metadata_noneg, permutations = 9999, by = "terms", na.action = na.omit)
#Df SumOfSqs      R2      F Pr(>F)
#SOM       1   0.2959 0.04511 0.8975  0.504
#Residual 19   6.2652 0.95489              
#Total    20   6.5612 1.00000           

# Same result.

# What if we include pH and GWC?

adonis2(ASV_tableRA.mat_o[1:24,] ~ pH_Ca + Moisture, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)

#Df SumOfSqs      R2      F Pr(>F)    
#pH_Ca     1   2.0743 0.27493 8.5617 0.0001 ***
#Moisture  1   0.3828 0.05074 1.5800 0.1127    
#Residual 21   5.0880 0.67434                  
#Total    23   7.5451 1.00000                  
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# so we should include only site in our permanova, or then pH without the site.
# why is this and we needed to check how the results change?

adonis2(ASV_tableRA.mat_o[1:24,] ~ pH_Ca , data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)

#           Df SumOfSqs      R2      F Pr(>F)    
#pH_Ca     1   2.0743 0.27493 8.3417  1e-04 ***
#Residual 22   5.4708 0.72507                  
#Total    23   7.5451 1.00000                  
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(ASV_tableRA.mat_o[1:24,] ~ Site, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)

#Df SumOfSqs      R2      F Pr(>F)    
#Site      3   4.2844 0.56784 8.7597  1e-04 ***
#Residual 20   3.2607 0.43216                  
#Total    23   7.5451 1.00000                  
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# which one we should have, Site or pH?
# If we have time, we could discuss "what Johanna and/or Antti would do in real life", i.e. if this would be an actual study.

#####################

# Influence of the site on alpha diversity:

# We can compare the diversities of the sites with for instance a t-test (remember the normality):

library(ggpubr)

ggplot(MMB117metadata_noneg, aes(x=Site, y=ASV_divSha)) +
  geom_point(aes(fill=factor(Site)), size=3, shape=21, colour="grey20",alpha=0.7,
             position=position_jitter(width=0.01, height=0.01)) +
  geom_boxplot (outlier.colour = NA, fill=NA, colour="grey20") +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  theme(panel.border = element_rect(colour = "black"))+
  #theme(axis.text.x=element_blank())+
  scale_fill_manual(values=c("greenyellow", "forestgreen", "darkslategray4","darkkhaki"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position="none")+
  geom_hline(yintercept = mean(na.omit(MMB117metadata_noneg$ASV_divSha)), linetype = 2)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Forest", hide.ns = TRUE) 

# so significantly higher in Gas station and park than in Forest.
# What does "significantly" mean?

# we can also get the comparison in a table format like this:

compare_means(ASV_divSha ~ Site,  data = MMB117metadata_noneg,
              ref.group = "Forest", method = "t.test")

# A tibble: 3 × 8
#.y.        group1 group2              p    p.adj p.format p.signif method
#<chr>      <chr>  <chr>           <dbl>    <dbl> <chr>    <chr>    <chr> 
#1 ASV_divSha Forest Field       0.0755    0.076    0.07552  ns       T-test
#2 ASV_divSha Forest Gas_station 0.0000273 0.000082 2.7e-05  ****     T-test
#3 ASV_divSha Forest Park        0.000154  0.00031  0.00015  ***      T-test


adonis2(ASV_tableRA.mat_o[1:24,] ~ Site + pH_Ca + SOM + Moisture, data=MMB117metadata_noneg, permutations=999, by = "terms",sqrt.dist=TRUE,na.action = na.omit)


ASVeucl<-dist(ASV_tableRA.mat_o[1:24,])

adonis2(ASVeucl ~ Site + pH_Ca + SOM + Moisture, data=MMB117metadata_noneg, permutations=999, by = "terms",na.action = na.omit)



MMB117anovaShan <- aov(MMB117metadata_noneg$ASV_divSha~ MMB117metadata_noneg$Site)
MMB117anovaShan

TukeyHSD(MMB117anovaShan)


ASVBray<-vegdist(ASV_tableRA.mat_o[1:24,])

length(ASVBray)

ASVBray

with(MMB117metadata_noneg,adonis2(ASVBray ~  pH_Ca + SOM + Moisture, data=MMB117metadata_noneg, permutations=999, by = "terms",na.action = na.omit))

rownames(ASV_tableRA.mat_o[1:24,])






#########################################

# So what you will write in your reports?



