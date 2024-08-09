##---------------------------------------------------------------
##Functions_MbThick.R
## Author: Doris Gomez
## Date Created: 2024-08-09
## Notes: this script lists all functions for the analyses presented in the article regarding membrane thickness
## Email: doris.gomez@cefe.cnrs.fr
##
##---------------------------------------------------------------




#######################FUNCTION TO SCALE VARIABLES EXCEPT MEMBRANE THICKNESS#######################################
ScaleContinuousVariables<-function(datatable,listvariablestoscale){
  datatable1<-datatable
  lengthlist<-length(as.vector(listvariablestoscale))
  for (j in 1:lengthlist){
    ifam<-match(listvariablestoscale[j],names(datatable))
    if (is.numeric(datatable[,ifam])){
      datatable1[,ifam]<-scale(datatable[,ifam])
    }
  }
  return(datatable1)

}



##################################EXPORTING OUTPUT STATISTICS###############################################

#function to fill the output table for three models
#adds the names of the factors and the name of the variable
FullTableFill<-function(A1,A2,A3,nbfacA1,nbfacA2,nbfacA3,typemodel="full"){
  # nbfacA1<-4; nbfacA2<-4; nbfacA3<-4
  # A1<-C1;A2<-C2;A3<-C3
  nbfactot<-nbfacA1+nbfacA2+nbfacA3
  tab <- matrix(character(0), nrow = nbfactot,ncol=7)
  for (j in 1:nbfacA1){
    tab[j,1]<-all.vars(formula(A1))[1]
    tab[j,2]<-rownames(A1$coefficients)[j]
    tab[j,3]<-paste0(round(A1$coefficients[j,1],4),"_",round(A1$coefficients[j,2],4))
    tab[j,4]<-round(A1$coefficients[j,3],3)
    tab[j,5]<-A1$coefficients[j,4]
    tab[j,6]<-A1$param["lambda"]
    tab[j,7]<-typemodel
  }
  for (j in 1:nbfacA2){
    tab[j+nbfacA1,1]<-all.vars(formula(A2))[1]
    tab[j+nbfacA1,2]<-rownames(A2$coefficients)[j]
    tab[j+nbfacA1,3]<-paste0(round(A2$coefficients[j,1],4),"_",round(A2$coefficients[j,2],4))
    tab[j+nbfacA1,4]<-round(A2$coefficients[j,3],3)
    tab[j+nbfacA1,5]<-A2$coefficients[j,4]
    tab[j+nbfacA1,6]<-A2$param["lambda"]
    tab[j+nbfacA1,7]<-typemodel
  }
  for (j in 1:nbfacA3){
    tab[j+nbfacA1+nbfacA2,1]<-all.vars(formula(A3))[1]
    tab[j+nbfacA1+nbfacA2,2]<-rownames(A3$coefficients)[j]
    tab[j+nbfacA1+nbfacA2,3]<-paste0(round(A3$coefficients[j,1],4),"_",round(A3$coefficients[j,2],4))
    tab[j+nbfacA1+nbfacA2,4]<-round(A3$coefficients[j,3],3)
    tab[j+nbfacA1+nbfacA2,5]<-A3$coefficients[j,4]
    tab[j+nbfacA1+nbfacA2,6]<-A3$param["lambda"]
    tab[j+nbfacA1+nbfacA2,7]<-typemodel
  }
  return (tab)
}

######################################## FUNCTIONS FOR GRAPHS   ###################################################################

#graph of the variation of the wing membrane thickness in the transparent zone MTT
GraphMTT<-function(mbg, intercept,slope_Transmittance,slope_FWLength,slope_PropTransp){
  # intercept<-A1$coefficients[1,1]
  # slope_Transmittance<-A1$coefficients[2,1]
  # slope_FWLength<-A1$coefficients[3,1]
  # slope_PropTransp<-A1$coefficients[4,1]
  median_Transmittance<-median(mbg$Transmittance)
  median_FWLength<-median(mbg$FWLength)
  median_PropTransp<-median(mbg$PropTransp)
  x1<-min(mbg$Transmittance)
  x2<-max(mbg$Transmittance)
  y1<-slope_Transmittance*x1+slope_FWLength*median_FWLength+slope_PropTransp*median_PropTransp+intercept
  y2<-slope_Transmittance*x2+slope_FWLength*median_FWLength+slope_PropTransp*median_PropTransp+intercept

  p<-ggplot(mbg, aes(x = Transmittance, y=MTT/1000,color=FWLength)) +
    geom_point(size=3) +theme_bw()+
    labs(y = "MTT = Mb thickness of the transparent zone (µm)",
         x = "Mean light transmittance over 300-700nm (%)")+
    labs(color='Forewing \nlength (mm)')+ylim(0,2.3)+
    theme(panel.grid = element_blank(),axis.text=element_text(size=14),
          axis.title=element_text(size=16),legend.title = element_text(size=14),legend.text = element_text(size=12),
          legend.position = "right",plot.title = element_text(size=18,hjust = 0.05))+
    geom_segment(x=x1,y=y1,xend=x2,yend=y2,lty=1,col="black",linewidth=1)

  return(p)}


#graph of the variation of the wing membrane thickness in the opaque zone MTO
GraphMTO<-function(mbg,intercept,slope_Transmittance,slope_FWLength,slope_PropTransp){
  # intercept<-A1$coefficients[1,1]
  # slope_Transmittance<-A1$coefficients[2,1]
  # slope_FWLength<-A1$coefficients[3,1]
  # slope_PropTransp<-A1$coefficients[4,1]
  median_Transmittance<-median(mbg$Transmittance)
  median_FWLength<-median(mbg$FWLength)
  median_PropTransp<-median(mbg$PropTransp)
  x1<-min(mbg$PropTransp)
  x2<-max(mbg$PropTransp)
  y1<-slope_Transmittance*median_Transmittance+slope_FWLength*median_FWLength+slope_PropTransp*x1+intercept
  y2<-slope_Transmittance*median_Transmittance+slope_FWLength*median_FWLength+slope_PropTransp*x2+intercept

  p<-ggplot(mbg, aes(x = PropTransp, y=MTO/1000,color=FWLength)) +
    geom_point(size=3) +theme_bw()+
    labs(y = "MTO = Mb thickness of the opaque zone (µm)",
         x = "Proportion of wing surface occupied by transparency (%)")+
    labs(color='Forewing \nlength (mm)')+ylim(0,3)+
    theme(panel.grid = element_blank(),axis.text=element_text(size=14),
          axis.title=element_text(size=16),legend.title = element_text(size=14),legend.text = element_text(size=12),
          legend.position = "right",plot.title = element_text(size=18,hjust = 0.05))+
    geom_segment(x=x1,y=y1,xend=x2,yend=y2,lty=1,col="black",linewidth=1)
  return(p)}


#graph of the variation of the difference in wing membrane thickness between the opaque zone and  the transparent zone DMT
#DMT=opaque - transparent
GraphDMT<-function(mbg,intercept,slope_Transmittance,slope_FWLength,slope_PropTransp){
  # intercept<-A1$coefficients[1,1]
  # slope_Transmittance<-A1$coefficients[2,1]
  # slope_FWLength<-A1$coefficients[3,1]
  # slope_PropTransp<-A1$coefficients[4,1]
  median_Transmittance<-median(mbg$Transmittance)
  median_FWLength<-median(mbg$FWLength)
  median_PropTransp<-median(mbg$PropTransp)
  x1<-min(mbg$PropTransp)
  x2<-max(mbg$PropTransp)
  y1<-slope_Transmittance*median_Transmittance+slope_FWLength*median_FWLength+slope_PropTransp*x1+intercept
  y2<-slope_Transmittance*median_Transmittance+slope_FWLength*median_FWLength+slope_PropTransp*x2+intercept

  p<-ggplot(mbg, aes(x = PropTransp, y=DMT/1000,color=FWLength)) +
    geom_point(size=3) +theme_bw()+
    labs(y = "DMT = Difference in Mb thickness (µm)",
         x = "Proportion of wing surface occupied by transparency (%)")+
    labs(color='Forewing \nlength (mm)')+ylim(-0.5,1.5)+
    theme(panel.grid = element_blank(),axis.text=element_text(size=14),
          axis.title=element_text(size=16),legend.title = element_text(size=14),legend.text = element_text(size=12),
          legend.position = "right",plot.title = element_text(size=18,hjust = 0.05))+
    geom_hline(aes(yintercept=0), colour="red", linetype="dashed",linewidth=1)+
    geom_segment(x=x1,y=y1,xend=x2,yend=y2,lty=1,col="black",linewidth=1)
  return(p)}
