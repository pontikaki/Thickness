##---------------------------------------------------------------
## make.R
## Author: Doris Gomez
## Date Created: 2024-08-09
## Notes: this script uses all functions and computes all steps of analyses and outputs presented in the article regarding membrane thickness
## Email: doris.gomez@cefe.cnrs.fr
##
##---------------------------------------------------------------



#### Clear memory ####
rm(list = ls())


########################################################## goes to description file to install and load all the packages needed
# downloads the external dependencies required by the project
devtools::install_deps()
#loads external dependencies and R functions
devtools::load_all()



#########################DATA READING#####################################################################
# reads phylogenetic tree
tree<-read.nexus(here::here("data","ithomiini_nexsansrien.tre"))
#reads input file with data on membrane thickness
mb<-data.frame(read.delim(here::here("data","TranspMbThickness.txt"),header=TRUE,stringsAsFactors = TRUE))


####################################################### LOADS THE FUNCTIONS NEEDED FOR THE ANALYSES
source(here::here("R","Functions_MbThick.R"))


####################################################### EXECUTION OF THE ANALYSES
source(here::here("analysis","Analyses_MbThick.R"))
