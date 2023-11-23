# script for the workshop on imputation of missing data in 
# functional traits databases
# mariasole calbi 1 dicember 2023, Pisa.

############1. SPECIES DATA / TRAITS GENERATION ##############
#install packages
install.packages("remotes")
library(remotes)
install.packages('traits') # run once
if(packageVersion("traits") == '0.2.0'){
devpars::install_github('terraref/traits', force = TRUE)}
remotes::install_github("ecoinfor/U.Taxonstand")
remotes::install_github("jinyizju/V.PhyloMaker2")
install.packages("http://cran.r-project.org/src/contrib/Archive/BHPMF/BHPMF_1.0.tar.gz", repos = NULL, type="source")

#load libraries
library(traits)
library(rgbif)
library(sf)
library(U.Taxonstand)
library(tidyr)
library(dplyr)
library(V.Phylomaker2)
library(ape)
library(stringr)
library(corrplot)
library(collapse)
library(mice)
library(VIM)
library(missForest)
library(PVR)
library(Rphylopars)
library(BHPMF)


setwd("~/Documents/PISA/script")#set wd

#get list of species from site of interest
#add site coordinates, orto botanico Pisa: 43.719988, 10.396481
PO_df <- data.frame(lon = c(10.396481),
                    lat = c(43.719988),
                    site= c("Pisa_orto"))

# Convert the dataframe to a spatial object. Note that the
# crs= 4326 parameter assigns a WGS84 coordinate system to the 
# spatial object
PO_sf <- st_as_sf(PO_df, coords = c("lon", "lat"), crs = 4326) 
cPO<-as.numeric(unlist(PO_sf$geometry))

Plantae <- 6#gbif key

dfcenter <- list("Pisa_orto"= cPO)#centroid of our site

Buffer <- 10000 #buffer in meters around site

My_Pol<- lapply(dfcenter,function(x){
  City_center <- st_sfc(st_point(x), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
  circle <- st_buffer(City_center, Buffer) 
  My_Polygon <- st_as_text(circle)
  ;My_Polygon}) #this function works also for multiple sites/centroids

gbif_taxon_keys <- c(6) #here one can also add animals if wanted
user="" # GBIF user name insert
pwd="" # GBIF password insert
email="" # your email insert
My_Pol
lapply(My_Pol[[1]],function(x){ #function to download gbif occurrences around centroid(s)
  gbif_download = occ_download(
    type="and",
    pred_in("taxonKey", gbif_taxon_keys),
    pred("hasGeospatialIssue", FALSE),
    pred("hasCoordinate", TRUE),
    pred_gte("year", 1980),
    pred_not(pred("basisOfRecord", "FOSSIL_SPECIMEN")),
    pred_within(x),
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email
  )
})

occ_download_wait('0056977-231002084531237') #substitute with download key obtained from previous function to get info on the status of the request
GBIF_PO <- occ_download_get('0056977-231002084531237') %>% #substitute with download key obtained from previous function to download the request
  occ_download_import()

#get trait data based on list of species
species_list<-unique(GBIF_PO$species)

#taxonomy check
#check nomenclature according to TPL
# load the TPL database
setwd("~/Documents/PISA/input/PHYLO_TPL")

library(openxlsx)
dat1 <- read.xlsx("Plants_TPL_database_part1.xlsx")
dat2 <- read.xlsx("Plants_TPL_database_part2.xlsx")
dat3 <- read.xlsx("Plants_TPL_database_part3.xlsx")
database <- rbind(dat1, dat2, dat3)
rm(dat1, dat2, dat3) #remove separate dfs

# run the function of name matching
res <- nameMatch(spList=species_list, spSource=database, author=FALSE, max.distance=1, Append=FALSE)

# save the results
setwd("~/Documents/PISA/input")
write.csv(res,"spec_ortoPisa_TPL.csv")
res<-read.csv("spec_ortoPisa_TPL.csv", header = T, row.names = 1)
#extract species with no Fuzzy matching
res2<-res[res$Fuzzy == FALSE,] #1792 spp

##LEDA traits retrieval

all <- c("canopy_height", "ldmc_geo", "leaf_mass", #list of traits of interest
         "leaf_size", "seed_mass",
         "seed_number")
out <- list() #create a list where to put selected traits data as elements
for (i in seq_along(all)) {
  cat(all[i], sep="\n")
  out[[i]] <- leda(all[i])
}
sapply(out, NROW) # number of trait entries in LEDA: 4934  3434  4470  4229 11952   687

sp_list_leda<-c(out[[1]]$sbs_name, out[[2]]$sbs_name, out[[3]]$sbs_name, 
                out[[4]]$sbs_name,out[[5]]$sbs_name, out[[6]]$sbs_name)
sp_list_leda<-unique(sp_list_leda) #count total spp covered: 3380

#taxonomy check for LEDA traits too
# run the function of name matching
resLEDA <- nameMatch(spList=sp_list_leda, spSource=database, author=FALSE, max.distance=1, Append=FALSE)
#extract species with no Fuzzy matching
resLEDA2<-resLEDA[resLEDA$Fuzzy == FALSE,]

sp_LEDA.oPisa<-res2[res2$Accepted_SPNAME %in% resLEDA2$Accepted_SPNAME,] #extract spp in common between our site data and LEDA data

#filter traits to this species pool
out.f<-lapply(out, function(df){
  df %>%
    filter(sbs_name %in% sp_LEDA.oPisa$Accepted_SPNAME)
})

canopy_height<-out.f[[1]][,c(1,4)] #shape traits data better to merge
colnames(canopy_height)<-c("sbs_name","canopy_height_m")
ldmc_geo<-out.f[[2]][,c(1,4)]
colnames(ldmc_geo)<-c("sbs_name","ldmc_geo_mgg")
leaf_mass<-out.f[[3]][,c(1,4)]
colnames(leaf_mass)<-c("sbs_name","leaf_mass_mg")
leaf_size<-out.f[[4]][,c(1,4)]
colnames(leaf_size)<-c("sbs_name","leaf_size_mm2")
seed_mass<-out.f[[5]][,c(1,5)]
colnames(seed_mass)<-c("sbs_name","seed_mass_mg")
seed_number<-out.f[[6]][,c(1,4)]
colnames(seed_number)<-c("sbs_name","seed_number")

traitsLEDA<-bind_rows(canopy_height, ldmc_geo, leaf_mass, leaf_size, #merge traits data
                  seed_mass,seed_number)

traitsLEDA.mean<-traitsLEDA %>% group_by(sbs_name) %>% summarize_at(c("canopy_height_m","ldmc_geo_mgg","leaf_mass_mg", "leaf_size_mm2", "seed_mass_mg",
                                    "seed_number"), mean, na.rm = TRUE) #get mean trait value by species since we have multiple obs by species in LEDA

traitsLEDA.mean[traitsLEDA.mean == "NaN"]<-NA  #use NA not NaN
options(scipen=999) #remove scientific notation it can cause problems later on
write.csv(traitsLEDA.mean, "LEDAtraits1dec.csv") #save file

#how many NAs are there in our LEDA-site data?
(sum(is.na(traitsLEDA.mean[,-1]))/prod(dim(traitsLEDA.mean[,-1])))*100 #34.38036
colSums(is.na(traitsLEDA.mean))
# sbs_name canopy_height_m    ldmc_geo_mgg    leaf_mass_mg   leaf_size_mm2    seed_mass_mg 
# 0             108             272             325             294             115 
# seed_number 
# 889       

LEDAnoNA<-na.omit(traitsLEDA.mean[,-7])#complete cases spp: 596 species for the 5 traits
sp_list<-as.vector(LEDAnoNA$sbs_name) #get complete cases spp list
LEDAnonar2<-res2[res2$Submitted_Name %in% LEDAnoNA$sbs_name,]

sp_list <- data.frame(species = LEDAnonar2$Accepted_SPNAME, genus = LEDAnonar2$Submitted_Genus, family = LEDAnonar2$Family)
sp_list <-sp_list[!duplicated(sp_list),]
LEDAnoNA<-LEDAnoNA[LEDAnoNA$sbs_name %in% sp_list$species,]

LEDAnonar2 <-LEDAnonar2[!duplicated(LEDAnonar2[,c(2:6)]),]
LEDAnoNA<-LEDAnoNA[LEDAnoNA$sbs_name %in% LEDAnonar2$Submitted_Name,]
sp_list$family[sp_list$family == "Compositae"]<-"Asteraceae" #change name to match megatree
sp_list$family[sp_list$family == "Leguminosae"]<-"Fabaceae" #change name

write.csv(LEDAnoNA, "LEDAcomplete_1dec.csv") #save
write.csv(LEDAnonar2, "LEDAcomplete_1dec_TAXOTPL.csv") #save

############2. PHYLOGENY ##############

# generate a phylogeny by pruning a mega tree for our species list
tree <- phylo.maker(sp.list = sp_list, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios = "S3")
write.tree(tree$scenario.3, "LEDAsptree.tre") #save

#visualize tree (circle form looks cooler)
plot(tree$scenario.3, type = "fan",
     label.offset = 1, cex = 0.1)

LEDA.ord<-as.data.frame(LEDAnoNA[order(match(LEDAnoNA$sbs_name, #reorder traits to match phylogeny tips order
                                               tree$tip.label)), ])
rownames(LEDA.ord)<-LEDA.ord$sbs_name
LEDA.ord$genus<-word(LEDA.ord$sbs_name, 1)#extract genus as first word of species name

#testing for traits normality 
shapiro.test(LEDA.ord$canopy_height_m)
shapiro.test(LEDA.ord$ldmc_geo_mgg)
shapiro.test(LEDA.ord$leaf_mass_mg)
shapiro.test(LEDA.ord$leaf_size_mm2)
shapiro.test(LEDA.ord$seed_mass_mg) #all non normally distribuited so a z and log transformations are needed

#test for correlation among traits using spearman because these are not normal
cort<-cor(as.matrix(LEDA.ord[,-c(1,7)]), method = "spearman") #non normal so use spearman
corrplot(cort) #some evident collinearities

#log and Z transform traits

#log 10 and z transform data
# logx <- log10(x)
# mlogx <- mean(logx, na.rm = T)
# slogx <- sd(logx, na.rm = T)
# x <- (logx - mlogx)/slogx # Z transformation

LEDA.ord.log<-log10(LEDA.ord[,-c(1,7)])#log 10 tranformed
LEDA.ord.logz<-as.data.frame(scale(LEDA.ord.log)) # log and z transformed will be used in the last imputation method
param<-scale(LEDA.ord.log)#keep track of back transformation parameters
attributes(param)#keep track of back transformation parameters 
LEDA.ord.logz$sbs_name<-LEDA.ord$sbs_name
LEDA.ord.logz$genus<-LEDA.ord$genus
write.csv(LEDA.ord.logz,"LEDA.ord.logz.csv") #save
write.csv(LEDA.ord, "LEDA.ord.csv") #save

#generate random NAs in variable percentages controlling to no introduce any row of only NAs
set.seed(123) #random seed very important to make analysis repeatable
LEDA10NA<-collapse::na_insert(LEDA.ord[,-c(1)], prop = .1) #insert 10% Nas randomly
LEDA10NA$genus<-LEDA.ord$genus
which(rowSums(is.na(LEDA10NA))==ncol(LEDA10NA)) #check if any row sums all NAs (we don't want that!)

set.seed(125)
LEDA30NA<-collapse::na_insert(LEDA.ord[,-c(1)], prop = .3)
LEDA30NA$genus<-LEDA.ord$genus
which(rowSums(is.na(LEDA30NA))==ncol(LEDA30NA))

write.csv(LEDA10NA, "LEDA10NA.csv") #save
write.csv(LEDA30NA, "LEDA30NA.csv") #save

#visualization of missing data
md.pattern(LEDA10NA)
md.pattern(LEDA30NA)

############3. IMPUTATION ##############
rm(list = ls()) #remove all elements from environment
options(scipen=999) #remove scientific notation
#load data
setwd("~/Documents/PISA/input")

LEDAnonar2<-read.csv("LEDAcomplete_1dec_TAXOTPL.csv", header = T, row.names = 1)
LEDA.ord.logz<-read.csv("LEDA.ord.logz.csv", header = T, row.names = 1)
LEDA.ord<-read.csv("LEDA.ord.csv", header = T, row.names = 1)#to be used in BHPMF
LEDA.ord<-LEDA.ord[,-1]
tree<-read.tree("LEDAsptree.tre")
LEDA10NA<-read.csv("LEDA10NA.csv", header = T, row.names = 1)
LEDA30NA<-read.csv("LEDA30NA.csv", header = T, row.names = 1)

setwd("~/Documents/PISA/output") #here we will store our outputs

#substitute missing values with column mean
NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)) #custom function for mean replacement of NAs
LEDA10NA.mean.all<-replace(LEDA10NA[,-6], TRUE, lapply(LEDA10NA[,-6], NA2mean))
# #or...
# LEDA10NA.mean.all<-as.data.frame(LEDA10NA) %>% 
#   mutate_all(funs(ifelse(is.na(.), mean(., na.rm = TRUE),.)))

LEDA30NA.mean.all<-replace(LEDA30NA[,-6], TRUE, lapply(LEDA30NA[,-6], NA2mean))#same in the 30% NAs dataset

#now we incorporate taxonomy as genera grouping
#substitute missing values with mean by group (genus), if not possible use general column mean
LEDA10NA.mean.g<-LEDA10NA %>%  #usinf dplyr and group_by
  group_by(genus) %>% 
  mutate_all(funs(ifelse(is.na(.), mean(., na.rm = TRUE),.)))%>% 
  ungroup()%>%
  mutate_all(funs(ifelse(is.na(.), mean(., na.rm = TRUE),.)))
LEDA10NA.mean.g<-as.data.frame(LEDA10NA.mean.g)

LEDA10NA.mean.g[LEDA10NA.mean.g == "NaN"]<-NA #change NaN for NA
rownames(LEDA10NA.mean.g)<-rownames(LEDA10NA)

#same for 30% NAs data
LEDA30NA.mean.g<-LEDA30NA %>% 
  group_by(genus) %>% 
  mutate_all(funs(ifelse(is.na(.), mean(., na.rm = TRUE),.))) %>% 
  ungroup()%>%
  mutate_all(funs(ifelse(is.na(.), mean(., na.rm = TRUE),.)))
LEDA30NA.mean.g<-as.data.frame(LEDA30NA.mean.g)

LEDA30NA.mean.g[LEDA30NA.mean.g == "NaN"]<-NA
rownames(LEDA30NA.mean.g)<-rownames(LEDA30NA)

#find missing data position in NAs dataframes to compare imputed and actual data
missing10<-which(is.na(LEDA10NA), arr.ind=TRUE) #position of NAs by rows and columns
missing30<-which(is.na(LEDA30NA), arr.ind=TRUE)

#Define the custom function that calculates RMSE (root mean squared error)
RMSE <- function(true_value, predicted_value){
  sqrt(mean((true_value - predicted_value)^2, na.rm = T))
}

#correlate values of imputed data with original complete cases data values
cor(as.numeric(LEDA.ord[missing10]), as.numeric(LEDA10NA.mean.all[missing10])) #0.3224331
#calculate RMSE of imputed data vs original complete cases data values
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(LEDA10NA.mean.all[missing10])) #9399.378

#same for 30% NAs dataset
cor(as.numeric(LEDA.ord[missing30]),as.numeric(LEDA30NA.mean.all[missing30]),use="complete.obs")#0.2869472
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(LEDA30NA.mean.all[missing30])) # 4699.22

#now for mean by genus groups + column mean
cor(as.numeric(LEDA.ord[missing10]), as.numeric(LEDA10NA.mean.g[missing10])) #0.3374773
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(LEDA10NA.mean.g[missing10])) # 9250.808

cor(as.numeric(LEDA.ord[missing30]), as.numeric(LEDA30NA.mean.g[missing30])) # 0.5869843
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(LEDA30NA.mean.g[missing30])) #3967.686

#kNN

imp_kNN_10<-kNN(LEDA10NA)  #default options, k = 5
imp_kNN_10<-imp_kNN_10[,1:5]
imp_kNN_30<-kNN(LEDA30NA)
imp_kNN_30<-imp_kNN_30[,1:5]

#check correlation and RMSE
cor(as.numeric(LEDA.ord[missing10]), as.numeric(imp_kNN_10[missing10])) # 0.859135
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(imp_kNN_10[missing10])) #7279.704

cor(as.numeric(LEDA.ord[missing30]), as.numeric(imp_kNN_30[missing30])) # 0.7412652
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(imp_kNN_30[missing30])) #3567.49

#Conditional Multiple Imputation: mice
options(scipen=999) #remove scientific notation
set.seed(123)
LD10.mice<-as.data.frame(complete(mice(LEDA10NA[,-6]))) #skips steps but can be unboxed for clarity of course; default options
LD30.mice<-as.data.frame(complete(mice(LEDA30NA[,-6])))

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing10]), as.numeric(LD10.mice[missing10]))#0.9429916
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(LD10.mice[missing10])) #3560.425

cor(as.numeric(LEDA.ord[missing30]),as.numeric(LD30.mice[missing30]))# 0.5156306
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(LD30.mice[missing30])) #4406.173

#random forest: missforest

set.seed(123)
missForest_impLD10 <-missForest(LEDA10NA[,-6])$ximp #default options

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing10]), as.numeric(missForest_impLD10[missing10]))#0.8798984
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(missForest_impLD10[missing10])) #5894.262

set.seed(123)#same for 30% NAs data
missForest_impLD30 <-missForest(LEDA30NA[,-6])$ximp

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing30]), as.numeric(missForest_impLD30[missing30]))#0.7625995
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(missForest_impLD30[missing30])) #3378.371

#####IMPUTATION WITH PHYLOGENETIC INFORMATION#####
#missForest with Phylogenetic information in the form of PVRs, Rphylopars, BHMPF

##calculate phylo signal of traits PVRs

LEDA10NA$tip<-gsub(" ", "_", rownames(LEDA10NA)) #obtain spp names that match tree tips format
LEDA30NA$tip<-gsub(" ", "_", rownames(LEDA30NA))
                   
PVRs<-PVRdecomp(tree, scale = TRUE) #PVR decomposition : Pcoa extracting eigenvectors
Eigen<-PVRs@Eigen #extrac eigenvectors
barplot(Eigen$values) #visualize eigen values

vectors<-as.data.frame(Eigen$vectors[,1:25])#keep first 25 eigvectors
tip<-PVRs@phylo$tip.label
PVRvect<-cbind(tip,vectors) #obtain PVRs and spp df to merge with traits data

#add PVR to traits
LEDA10NA_PVR<-merge(LEDA10NA , PVRvect, by="tip",all.x=TRUE) #merge traits data with PVRs
LEDA10NA_PVR<-LEDA10NA_PVR[,-7]

#missForest workflow
set.seed(123)
traits_f_PVR.imp10 <- missForest(LEDA10NA_PVR[,-1], verbose = TRUE)$ximp #imputation

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing10]), as.numeric(traits_f_PVR.imp10[missing10]))# 0.8050689
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(traits_f_PVR.imp10[missing10])) # 6717.548

#now same with 30% NAs data
#add PVR to traits
LEDA30NA_PVR<-merge(LEDA30NA , PVRvect, by="tip",all.x=TRUE) #merge PVRs with traits data
LEDA30NA_PVR<-LEDA30NA_PVR[,-7]

traits_f_PVR.imp30 <- missForest(LEDA30NA_PVR[,-1], verbose = TRUE)$ximp #imputation

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing30]), as.numeric(traits_f_PVR.imp30[missing30]))# 0.5585047
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(traits_f_PVR.imp30[missing30])) # 4062.18

## RPhylopars

LEDA10NA$species<-LEDA10NA$tip# prep df
LEDA10NA<-LEDA10NA[,-c(6,7)]
LEDA10NA<-LEDA10NA[,c(6,1,2,3,4,5)]

phypars10 <- phylopars(trait_data = LEDA10NA, tree = tree) #imputation with default options: Brownian motion
phypars10I<-as.data.frame(phypars10$anc_recon[1:549,])
rownames(LEDA10NA)<-LEDA10NA$species
phypars10I<-phypars10I[match(rownames(LEDA10NA), rownames(phypars10I)), ]

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing10]), as.numeric(phypars10I[missing10]))# 0.9498553
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(phypars10I[missing10]))# 3753.587

#same for 30% NAs data
LEDA30NA$species<-LEDA30NA$tip
LEDA30NA<-LEDA30NA[,-c(6,7)]
LEDA30NA<-LEDA30NA[,c(6,1,2,3,4,5)]

phypars30 <- phylopars(trait_data = LEDA30NA, tree = tree)
phypars30I<-as.data.frame(phypars30$anc_recon[1:549,])
rownames(LEDA30NA)<-LEDA30NA$species
phypars30I<-phypars30I[match(rownames(LEDA30NA), rownames(phypars30I)), ]

#also scatter plot is possibile to visualize correlation
plot(LEDA.ord[missing30],phypars30I[missing30],pch=21,bg="grey",cex=1.5,
     xlab="true values",ylab="imputed values 30 % NAs using Rphylopars",asp=1)
lines(x=c(min(par()$usr[c(1,3)]),max(par()$usr[c(2,4)])),
      y=c(min(par()$usr[c(1,3)]),max(par()$usr[c(2,4)])),
      col="blue")
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing30]), as.numeric(phypars30I[missing30]))#0.8890581
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(phypars30I[missing30]))# 2348.346

#try with different evolution models....
phypars10OU <- phylopars(trait_data = LEDA10NA, tree = tree, model = "OU") #with Ornstein-Uhlenbeck
phypars10OUI<-as.data.frame(phypars10OU$anc_recon[1:549,])
phypars10OUI<-phypars10OUI[match(rownames(LEDA10NA), rownames(phypars10OUI)), ]

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing10]), as.numeric(phypars10OUI[missing10]))#  0.9500217
RMSE(as.numeric(LEDA.ord[missing10]), as.numeric(phypars10OUI[missing10]))# 3771.924

#now with 30% NAs
phypars30OU <- phylopars(trait_data = LEDA30NA, tree = tree, model = "OU")
phypars30OUI<-as.data.frame(phypars30OU$anc_recon[1:549,])
phypars30OUI<-phypars30OUI[match(rownames(LEDA30NA), rownames(phypars30OUI)), ]

#check cor and RMSE
cor(as.numeric(LEDA.ord[missing30]), as.numeric(phypars30OUI[missing30]))# 0.874755
RMSE(as.numeric(LEDA.ord[missing30]), as.numeric(phypars30OUI[missing30]))# 2438.515

##bayesian BHMPF

#10%NA
hierarchy.info<-LEDAnonar2[,c(21,8,19)] #hierarchy df: sp, genus, family infos
colnames(hierarchy.info)<-c("species", "genus", "family")
hierarchy.info$species<-gsub(" ", "_", hierarchy.info$species) 
#order species according LEDA30NA
hierarchy.info<-hierarchy.info[order(match(hierarchy.info$species, LEDA10NA$species)),]
hierarchy.info<-hierarchy.info %>% mutate(plant_id = 1:n()) %>% select(plant_id, everything())

trait.info<-LEDA10NA[rowSums(is.na(LEDA10NA)) != ncol(LEDA10NA), ]#trait info should be a matrix in the same order as hierarchy.info
trait.info<-as.matrix(trait.info[,-1])

#log- and z-transformations and back to improve normality and equalize measurements in 
#the cost function during optimization
back_trans_pars <- list()
rm_col <- c()
for(i in 1:ncol(trait.info)){
  x <- trait.info[,i] # goes through the columns
  min_x <- min(x,na.rm = T) # takes the min of each column
  if(min_x < 0.00000000001){
    x <- x - min_x + 1 # make this optional if min x is neg
  }
  logx <- log10(x)
  mlogx <- mean(logx, na.rm = T)
  slogx <- sd(logx, na.rm = T)
  x <- (logx - mlogx)/slogx # Z transformation
  back_trans_pars[[i]] <- list(min_x = min_x,
                               mlogx = mlogx,
                               slogx = slogx)
  trait.info[,i] <- x
}
write.csv(back_trans_pars, "back_trans_pars.csv")#save

#gap filling with no hierarchy information
GapFilling(trait.info, hierarchy.info,
           used.num.hierarchy.levels = 0, tuning = T, 
           mean.gap.filled.output.path = "LEDAmean_gap_filled10_NOHYE.txt",
           std.gap.filled.output.path= "LEDAstd_gap_filled10_NOHYE.txt")

#gap filling with hierarchy information
GapFilling(trait.info, hierarchy.info, tuning = T,
           mean.gap.filled.output.path = "LEDAmean_gap_filled10.txt",
           std.gap.filled.output.path= "LEDAstd_gap_filled10.txt")

B.LEDA10NA_GFNOHYE<-read.table("LEDAmean_gap_filled10_NOHYE.txt", header = T)
rownames(B.LEDA10NA_GFNOHYE)<-rownames(LEDA10NA)

B.LEDA10NA_GF<-read.table("LEDAmean_gap_filled10.txt", header = T)
rownames(B.LEDA10NA_GF)<-rownames(LEDA10NA)

#check imputed values against full log/z-transformed traits values
cor(as.numeric(LEDA.ord.logz[missing10]), as.numeric(B.LEDA10NA_GFNOHYE[missing10]))#  0.1850964
RMSE(as.numeric(LEDA.ord.logz[missing10]), as.numeric(B.LEDA10NA_GFNOHYE[missing10]))#  1.566254

#for hierarchy info included gap filling
cor(as.numeric(LEDA.ord.logz[missing10]), as.numeric(B.LEDA10NA_GF[missing10]))#  0.5519345
RMSE(as.numeric(LEDA.ord.logz[missing10]), as.numeric(B.LEDA10NA_GF[missing10]))# 1.018367

#30%NA
trait.info30<-LEDA30NA[rowSums(is.na(LEDA30NA)) != ncol(LEDA30NA), ]
trait.info30<-as.matrix(trait.info30[,-1])

#log- and z-transformations and back to improve normality and equalize measurements in 
#the cost function during optimization
back_trans_pars30 <- list()
rm_col <- c()
for(i in 1:ncol(trait.info30)){
  x <- trait.info30[,i] # goes through the columns
  min_x <- min(x,na.rm = T) # takes the min of each column
  if(min_x < 0.00000000001){
    x <- x - min_x + 1 # make this optional if min x is neg
  }
  logx <- log10(x)
  mlogx <- mean(logx, na.rm = T)
  slogx <- sd(logx, na.rm = T)
  x <- (logx - mlogx)/slogx # Z transformation
  back_trans_pars30[[i]] <- list(min_x = min_x,
                               mlogx = mlogx,
                               slogx = slogx)
  trait.info30[,i] <- x
}
write.csv(back_trans_pars30, "back_trans_pars30.csv")

#gap filling no hierarchy information
GapFilling(trait.info30, hierarchy.info,
           used.num.hierarchy.levels = 0, tuning = T, 
           mean.gap.filled.output.path = "LEDAmean_gap_filled30_NOHYE.txt",
           std.gap.filled.output.path= "LEDAstd_gap_filled30_NOHYE.txt")

#gap filling with hierarchy information
GapFilling(trait.info30, hierarchy.info, tuning = T,
           mean.gap.filled.output.path = "LEDAmean_gap_filled30.txt",
           std.gap.filled.output.path= "LEDAstd_gap_filled30.txt")

B.LEDA30NA_GFNOHYE<-read.table("LEDAmean_gap_filled30_NOHYE.txt", header = T)
rownames(B.LEDA30NA_GFNOHYE)<-rownames(LEDA30NA)

B.LEDA30NA_GF<-read.table("LEDAmean_gap_filled30.txt", header = T)
rownames(B.LEDA30NA_GF)<-rownames(LEDA30NA)

#check imputed values against full logztransformed traits values
cor(as.numeric(LEDA.ord.logz[missing30]), as.numeric(B.LEDA30NA_GFNOHYE[missing30]))# 0.04600237
RMSE(as.numeric(LEDA.ord.logz[missing30]), as.numeric(B.LEDA30NA_GFNOHYE[missing30]))#  2.095398

cor(as.numeric(LEDA.ord.logz[missing30]), as.numeric(B.LEDA30NA_GF[missing30]))# 0.3074057
RMSE(as.numeric(LEDA.ord.logz[missing30]), as.numeric(B.LEDA30NA_GF[missing30]))# 1.471656

###GRAPHS of all methods correlations and RMSE####
results<-read.csv("cor_RMSEresults.csv", header = T, sep = ";", check.names = F)
barplot(results[results$`Nas_%` == "10",]$cor , main="correlations 10% NAs",ylab="cor",las=2, names.arg = results[results$`Nas_%` == "10",]$Method)
barplot(results[results$`Nas_%` == "30",]$cor , main="correlations 30% NAs",ylab="cor",las=2, names.arg = results[results$`Nas_%` == "30",]$Method )


