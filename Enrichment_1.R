#--------------------------------gprofiler2--------------------------------------------
#read data 
dat <- read.csv(file.choose())
#splite the accessions  as a list 
Acc  <- split(dat , dat$Accessions) 
#load gprofiler2 package 
install.packages("gprofiler2")
library("gprofiler2")
#used for better performance  
library("dplyr")
# use gost object to extract the data by Applying a Query of Accession number list 
#i used just 5 accessions due to runtime 
gostres <-gost(query = Acc[1:5],organism = "hsapiens")
# convert gostres list to a dataframe then select some columns   
gostres_df <- data.frame(gostres$result) %>% select("query","p_value"
                                                   ,"term_id","term_name",
                                                   "source","intersection_size")
# the same result as abve but different way without dplyr package
selected_data <- data.frame(cbind(gostres_df$query,gostres_df$p_value
                                  ,gostres_df$term_id,gostres_df$term_name,
                                  gostres_df$source,gostres_df$intersection_size))
# used for more description of columns 
colnames(selected_data)<- c("Accession","p_value","term_id","term_name","DataBase","intersection_size")



#-------------------------------------------------------------------------------
#--------------------------------String--------------------------------------------
#install the package 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STRINGdb")
#load the package 
library("STRINGdb")
#for more details read the documentation 
browseVignettes("STRINGdb")
# make an object called string_dp while have the data samples of the database
#functional network for loading full functions 
# or 
# physical for physical subnetwork which link only the proteins which share a physical complex 
string_db <- STRINGdb$new( version="11.5", species=9606,
      score_threshold=200, network_type="functional", input_directory=getwd())

#map the data to the object for adding more detials like string_id 
example1_mapped <- string_db$map( dat, "Accessions", removeUnmappedRows = TRUE )
#using string_id for plot the network for just first 100 protien 
hits <- example1_mapped$STRING_id[1:100]
# for show the network 
string_db$plot_network( hits )
#get enrichment fir the first 100 accessions 
Enrichment <- string_db$get_enrichment(hits)
#show some data 
head(Enrichment)
# show annotations 
annotations <- string_db$get_annotations( hits )
head(annotations, n=20)


#interaction between protiens (from) (to) (combined score) based on string_id 
interacted_protiens <-string_db$get_interactions(hits)

# get clusters
clustersList <- string_db$get_clusters(example1_mapped$STRING_id[1:100])
# plot first 4 clusters
par(mfrow=c(1,2))
for(i in seq(1:2)){
   string_db$plot_network(clustersList[[i]])
   }
#-------------------------------------------------------------------------------------------------------
#--------------------------------Enrichr--------------------------------------------
# transform protiens id to gene sympoles
library("UniprotR")
Accessions <-GetAccessionList("D:/Work57/NSAF.csv") 
#Get Taxanomy Information 
TaxaObj <- GetNamesTaxa(Accessions) 
#install erichR package  
install.packages("enrichR") 
#load package 
library("enrichR")
#list of ErivhR sites for all available organisms 
listEnrichrSites() 
#select just human enrichr 
setEnrichrSite("Enrichr")
# connect to website live  
websiteLive <- TRUE
#load all Enrichr databases 
#the list of all available databases from Enrichr
dbs <- listEnrichrDbs() 
#error handling: in case no databases found then cut the connection 
if (is.null(dbs)) websiteLive <- FALSE
# if there is a connection then show the head of databases  
if (websiteLive) head(dbs)
# choose selectecd databases 
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
# use gene names from taxa object from uniprotR 
if (websiteLive) {
  enriched <- enrichr(TaxaObj$Gene.Names, dbs)
}
#gat all the data about Biological process by genenames  
enriched_BP<- if (websiteLive) enriched[["GO_Biological_Process_2015"]]
# plot the Enrichment 
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
#---------------------------------------------------------------------------------------
#--------------------------------Go ontolgy--------------------------------------------

#Get Gene ontolgy Information 
GeneOntologyObj <- GetProteinGOInfo(Accessions) 
#Plot Biological process information top 10 go terms  
PlotGOBiological(GeneOntologyObj, Top = 10) 
#Plot molecular function information top 20 go terms
Plot.GOMolecular(GeneOntologyObj, Top = 20)
#Plot subcellualr localization information 
Plot.GOSubCellular(GeneOntologyObj) 
#Combine Gene ontology plots into one plot 
PlotGoInfo(GeneOntologyObj)
#Handy visualization for publications 
PlotGOAll(GOObj = GeneOntologyObj, Top = 10, directorypath = getwd(), width = 8, height = 5)





