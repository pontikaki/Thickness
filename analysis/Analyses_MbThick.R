##---------------------------------------------------------------
##Analyses_MbThick.R
## Author: Doris Gomez
## Date Created: 2024-08-09
## Notes: this script runs the analyses presented in the article regarding membrane thickness
## Email: doris.gomez@cefe.cnrs.fr
##
##---------------------------------------------------------------


#source(here::here("R","Functions_MbThick.R"))

names(mb)
list<-names(mb[,3:13])
list

mbs<-ScaleContinuousVariables(mb,list)
summary(mbs)


#######################PGLS analyses#######################################################################
# selection of the species with both transparent and opaque zones
mbsb<-mbs[(!is.na(mbs$FWLength))&(mbs$SpType=="B"),]
#
table(mbsb$DMTSign)
droplevels(mbsb$phylotag)

#we define first the tree on which to work
abc<-mbsb$phylotag
def<-factor(tree$tip.label)
dsq<-intersect(abc,def)
parttree1<-keep.tip(tree,dsq)

#correlation matrix between predictors
cor_matrix <- cor(mbsb[, c("Transmittance", "FWLength", "PropTransp")])
cor.test(mbsb$Transmittance,mbsb$FWLength)
cor.test(mbsb$PropTransp,mbsb$FWLength)
cor.test(mbsb$Transmittance,mbsb$PropTransp)
print(cor_matrix)

#correlation between the predictors, values between 1 and 2 are moderate correlations
lm_model <- lm(MTT/1000 ~ Transmittance+FWLength+PropTransp, data = mbsb)
vif_values <- car::vif(lm_model)
print(vif_values)


#define data table for pgls
data.comp1 <- caper::comparative.data(phy = parttree1, data = mbsb, vcv = TRUE, names.col = phylotag, na.omit = FALSE)

#data on membrane thickness MTT MTO DMT are expressed in nm but too large numbers,
#better to express them in micrometers, so divide by 1000 all the analyses and graphs

#best model for MTT =model 2
A<-list()
A[[1]]<-caper::pgls(formula = MTT/1000 ~  Transmittance+FWLength+PropTransp, data=data.comp1, lambda ="ML")
A[[2]]<-caper::pgls(formula = MTT/1000 ~  Transmittance+FWLength, data=data.comp1, lambda ="ML")
A[[3]]<-caper::pgls(formula = MTT/1000 ~  Transmittance+PropTransp, data=data.comp1, lambda ="ML")
A[[4]]<-caper::pgls(formula = MTT/1000 ~  FWLength+PropTransp, data=data.comp1, lambda ="ML")
A[[5]]<-caper::pgls(formula = MTT/1000 ~  Transmittance, data=data.comp1, lambda ="ML")
A[[6]]<-caper::pgls(formula = MTT/1000 ~  FWLength, data=data.comp1, lambda ="ML")
A[[7]]<-caper::pgls(formula = MTT/1000 ~  PropTransp, data=data.comp1, lambda ="ML")
A[[8]]<-caper::pgls(formula = MTT/1000 ~  1, data=data.comp1, lambda ="ML")
A[[1]]$aicc
A[[2]]$aicc
A[[3]]$aicc
A[[4]]$aicc
A[[5]]$aicc
A[[6]]$aicc
A[[7]]$aicc
A[[8]]$aicc

#bestmodel for MTO=model 4
B<-list()
B[[1]]<-caper::pgls(formula = MTO/1000 ~  Transmittance+FWLength+PropTransp, data=data.comp1, lambda ="ML")
B[[2]]<-caper::pgls(formula = MTO/1000 ~  Transmittance+FWLength, data=data.comp1, lambda ="ML")
B[[3]]<-caper::pgls(formula = MTO/1000 ~  Transmittance+PropTransp, data=data.comp1, lambda ="ML")
B[[4]]<-caper::pgls(formula = MTO/1000 ~  FWLength+PropTransp, data=data.comp1, lambda ="ML")
B[[5]]<-caper::pgls(formula = MTO/1000 ~  Transmittance, data=data.comp1, lambda ="ML")
B[[6]]<-caper::pgls(formula = MTO/1000 ~  FWLength, data=data.comp1, lambda ="ML")
B[[7]]<-caper::pgls(formula = MTO/1000 ~  PropTransp, data=data.comp1, lambda ="ML")
B[[8]]<-caper::pgls(formula = MTO/1000 ~  1, data=data.comp1, lambda ="ML")
B[[1]]$aicc
B[[2]]$aicc
B[[3]]$aicc
B[[4]]$aicc
B[[5]]$aicc
B[[6]]$aicc
B[[7]]$aicc
B[[8]]$aicc

#best model for DMT=model 7
D<-list()
D[[1]]<-caper::pgls(formula = DMT/1000 ~  Transmittance+FWLength+PropTransp, data=data.comp1, lambda ="ML")
D[[2]]<-caper::pgls(formula = DMT/1000 ~  Transmittance+FWLength, data=data.comp1, lambda ="ML")
D[[3]]<-caper::pgls(formula = DMT/1000 ~  Transmittance+PropTransp, data=data.comp1, lambda ="ML")
D[[4]]<-caper::pgls(formula = DMT/1000 ~  FWLength+PropTransp, data=data.comp1, lambda ="ML")
D[[5]]<-caper::pgls(formula = DMT/1000 ~  Transmittance, data=data.comp1, lambda ="ML")
D[[6]]<-caper::pgls(formula = DMT/1000 ~  FWLength, data=data.comp1, lambda ="ML")
D[[7]]<-caper::pgls(formula = DMT/1000 ~  PropTransp, data=data.comp1, lambda ="ML")
D[[8]]<-caper::pgls(formula = DMT/1000 ~  1, data=data.comp1, lambda ="ML")
D[[1]]$aicc
D[[2]]$aicc
D[[3]]$aicc
D[[4]]$aicc
D[[5]]$aicc
D[[6]]$aicc
D[[7]]$aicc
D[[8]]$aicc

#we take the same full model to fill the table
C<-list()
C[[1]]<-caper::pgls(formula = MTT/1000 ~  Transmittance+FWLength+PropTransp, data=data.comp1, lambda ="ML")
C1<-summary(C[[1]])
C[[2]]<-caper::pgls(formula = MTO/1000 ~  Transmittance+FWLength+PropTransp, data=data.comp1, lambda ="ML")
C2<-summary(C[[2]])
C[[3]]<-caper::pgls(formula = DMT/1000 ~  Transmittance+FWLength+PropTransp, data=data.comp1, lambda ="ML")
C3<-summary(C[[3]])
# C1
# C2
# C3

#we take the same full model to fill the table
Q<-list()
Q[[1]]<-caper::pgls(formula = MTT/1000 ~  Transmittance+FWLength, data=data.comp1, lambda ="ML")
Q1<-summary(Q[[1]])
Q[[2]]<-caper::pgls(formula = MTO/1000 ~  FWLength+PropTransp, data=data.comp1, lambda ="ML")
Q2<-summary(Q[[2]])
Q[[3]]<-caper::pgls(formula = DMT/1000 ~  PropTransp, data=data.comp1, lambda ="ML")
Q3<-summary(Q[[3]])
# Q1
# Q2
# Q3



nbrowsfull<-4
#runs the function for the full models
tabfull<-FullTableFill(C1,C2,C3,nbrowsfull,nbrowsfull,nbrowsfull,typemodel="full")
#runs the function for the best models
tabbest<-FullTableFill(Q1,Q2,Q3,3,3,2,typemodel="best")
alpha<-as.data.frame(rbind(tabfull,tabbest))
names(alpha)<-c("dependent variable","factor","estimate_se","t-value","p-value","lambda","typemodel")

#writes the output statistical table
write.table(alpha,here::here("outputs","table1.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)





######################################## FUNCTIONS FOR GRAPHS   ###################################################################
#we select the table on which to plot the results, variables are not scaled
mbg<-mb[(!is.na(mb$FWLength))&(mb$SpType=="B"),]



################################# RUNS GRAPHS AND SAVES IN OUTPUT FOLDER #########################################
#intercept, slopeTransmittance, slopeLength, slopePropTrans

svg(file=here::here("outputs","Graph_MTTfull.svg"), width=10, height=5)
p<-GraphMTT(mbg, C1$coefficients[1,1],C1$coefficients[2,1],C1$coefficients[3,1],C1$coefficients[4,1])
print(p)
dev.off()
svg(file=here::here("outputs","Graph_MTOfull.svg"), width=10, height=5)
p<-GraphMTO(mbg, C2$coefficients[1,1],C2$coefficients[2,1],C2$coefficients[3,1],C2$coefficients[4,1])
print(p)
dev.off()
svg(file=here::here("outputs","Graph_DMTfull.svg"), width=10, height=5)
p<-GraphDMT(mbg, C3$coefficients[1,1],C3$coefficients[2,1],C3$coefficients[3,1],C3$coefficients[4,1])
print(p)
dev.off()


svg(file=here::here("outputs","Graph_MTTbest.svg"), width=10, height=5)
p<-GraphMTT(mbg, Q1$coefficients[1,1],Q1$coefficients[2,1],Q1$coefficients[3,1],0)
print(p)
dev.off()
svg(file=here::here("outputs","Graph_MTObest.svg"), width=10, height=5)
p<-GraphMTO(mbg, Q2$coefficients[1,1],0,Q2$coefficients[2,1],Q2$coefficients[3,1])
print(p)
dev.off()
svg(file=here::here("outputs","Graph_DMTbest.svg"), width=10, height=5)
p<-GraphDMT(mbg, Q3$coefficients[1,1],0,0,Q3$coefficients[2,1])
print(p)
dev.off()


