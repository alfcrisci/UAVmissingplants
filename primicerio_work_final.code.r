#############################################################
# Setup working directory in your desktop

setwd("/home/alf/Scrivania/lav_primicerio_final")

#############################################################
# Load dependencies and functions

# Before load library check and install the packages as usual if necessary.

#############################################################
# A Load spatial libraries

library(maptools)
library(spatstat)
library(spdep)
library(rgdal)
library(raster)
library(rgeos)
library(leaflet)
library(Cairo)
library(MASS)
library(rpart.plot)

citation("maptools")
citation("spatstat")
citation("spdep")
citation("rgdal")
citation("rpart.plot")
citation("spdep")
citation("rgeos")
citation("MASS")
citation("raster")

#############################################################
# B Load graphical libraries

library(ggplot2)
library(sjPlot)
library(sjmisc)

citation("ggplot2")
citation("sjPlot")

#############################################################
# C Load modeling libraries

library(rpart)
library(plotROC)


citation("rpart")
citation("plotROC")


#############################################################
# load other function

source("auxillary_bulichella.r")

################################################################################################
# Load geo data

mat_bulichella_sp=readRDS("geo_bulichella/mat_bulichella.rds")

buliGEO=brick("raster/buliGEO.tif")
buli_mask=brick("raster/C_mask.tif")

#########################################################################################################################à
# trattamenti raster per calcolo bound e distanze

r_bulichella <- stack(SpatialPixelsDataFrame(mat_bulichella_sp, tolerance = 0.00001, mat_bulichella_sp@data))
proj4string(r_bulichella)=CRS("+init=epsg:4326")
writeRaster(r_bulichella,"raster/r_bulichella.tif",overwrite=TRUE)

r_mask_bulichella=r_bulichella[["Area"]]*1
setValues(r_mask_bulichella)=NA
writeRaster(r_mask_bulichella,"raster/r_mask_bulichella.tif",overwrite=TRUE)


# Calculate Dist_bound parameter in unit images 
# IT Calcolo della forma convessa e della distanza

qhull_bulichella=gConvexHull(mat_bulichella_sp)
proj4string(qhull_bulichella)=CRS("+init=epsg:4326")
coords_edges_bulichella=as.data.frame(getEdges(qhull_bulichella))
coords_edges_bulichella$id=1
coordinates(coords_edges_bulichella)= ~ x+y
class(coords_edges_bulichella)
bound_line_bulichella=gBoundary(qhull_bulichella)
saveRDS(bound_line_bulichella,"geo_bulichella/bound_line_bulichella.rds")



bound_line_bulichella=readRDS("geo_bulichella/bound_line_bulichella.rds")
proj4string(bound_line_bulichella)=proj4string(mat_bulichella_sp)
dist_bound=gDistance(mat_bulichella_sp,bound_line_bulichella,byid=T) 


# Calculate factor variable ( 0- 1, False True ) quantile UnderPerc_Area if feature is under the 90th percentile 

q_area=quantile(values(r_bulichella[["Area"]]),probs = seq(0, 1, 0.1),na.rm=T)
mat_bulichella_sp$Underperc=extract(r_bulichella[["Area"]]>q_area[10],mat_bulichella_sp)


####################################################################################
# estract missing plants fetatures

mat_bulichella_miss=mat_bulichella_sp[mat_bulichella_sp$MISSING>0,]

##################################################################################
# Create a spatastat spatial object to visualize missing plant density


sSp_bulichella_miss <- as(SpatialPoints(mat_bulichella_miss), "ppp")  # convert points to pp class
Dens_bulichella_miss <- density(sSp_bulichella_miss, adjust = 0.2)  # create density object

class(Dens_bulichella_miss)  # just for interest: it's got it's of pixel image class


plot(Dens_bulichella_miss)  # default plot for spatial density
contour(density(sSp_bulichella_miss, adjust = 0.2), nlevels = 4)  # plot as contours - this is where we're heading plot of chunk Contour plot 

Dsg_bulichella_miss <- as(Dens_bulichella_miss, "SpatialGridDataFrame")  # convert to spatial grid class
Dim_bulichella_miss <- as.image.SpatialGridDataFrame(Dsg_bulichella_miss)  # convert again to an image
Dcl_bulichella_miss <- contourLines(Dim_bulichella_miss, nlevels = 8)  # create contour object - change 8 for more/fewer levels
SLDF_bulichella_miss<- ContourLines2SLDF(Dcl_bulichella_miss, CRS("+init=epsg:4326"))  # convert to SpatialLinesDataFrame
SLDF_bulichella_miss=SLDF_bulichella_miss[SLDF_bulichella_miss$level!=0,] # leave data boudary
plot(SLDF_bulichella_miss, col = terrain.colors(4))

################################################################################################################

CairoPNG(filename = "results/plot_bulichella_image.png",res=300)

plotRGB(buliGEO)
plot(bound_line_bulichella,col ='red',add=TRUE)

dev.off()

CairoPNG(filename = "results/plot_bulichella_over_mask.png",res=300)
plotRGB(buliGEO)
plotRGB(buli_mask,alpha=120,colNA='red',add=TRUE)
plot(bound_line_bulichella,col ='red',add=TRUE)
dev.off()


CairoPNG(filename = "results/plot_density_missing.png",res=300)
plotRGB(buliGEO)
plot(bound_line_bulichella,col ='red',add=TRUE)
plot(SLDF_bulichella_miss, col = terrain.colors(4),add=TRUE)
plot(mat_bulichella_miss,pch = 19,cex = .1,col ='brown2',add=TRUE)
dev.off()

saveRDS(SLDF_bulichella_miss,"geo_bulichella/SLDF_bulichella_miss.rds")


####################################################################################
# Local moran calculation by using of N matrix of neighbours


mat_bulichella_50=nb2listw(knn2nb(knearneigh(mat_bulichella_sp,k=50))) # 50 plants
mat_bulichella_30=nb2listw(knn2nb(knearneigh(mat_bulichella_sp,k=30))) # 30 plants


mat_bulichella_sp@data$L_moran_50=localmoran(mat_bulichella_sp$Area, mat_bulichella_50)[,1]
mat_bulichella_sp@data$L_moran_30=localmoran(mat_bulichella_sp$Area, mat_bulichella_30)[,1]
mat_bulichella_sp@data$L_moran_30_p=localmoran(mat_bulichella_sp$Perimetro, mat_bulichella_30)[,1]

###############################################################################################
# LOESS deviation by line vigour model : model_primicerio function available in auxillary_bulichella.r

mat_bulichella_ls=split(mat_bulichella_sp@data,mat_bulichella_sp@data$FILARE)

res=list()
for (i in seq_along(mat_bulichella_ls)) {
  res[[i]]=model_primicerio(mat_bulichella_ls[[i]],
                            saveplot=TRUE,
                            titlefig=paste0("Modeling Plant Missing FILARE_",as.character(i)),
                            namefig=paste0("results/Modeling_Plant_Missing_FILARE ",as.character(i),".png"),
                            treshshold=100)
}


ls_model_residuals=lapply(res,function(x) x$model_residuals)

mat_bulichella_sp@data$Line_res=unlist(ls_model_residuals)

#########################################################################################################
# create guess variable candidate of missing when deviation are higher than treshsold

ls_canditate=lapply(res,function(x) x$vector)

mat_bulichella_sp@data$candidate=unlist(ls_canditate)

names(mat_bulichella_sp@data)

saveRDS(mat_bulichella_sp,"mat_bulichella_sp.rds")


#############################################################################################################################
# Modeling steps

mat_bulichella_sp=readRDS("mat_bulichella_sp.rds")

#############################################################################################################################
# model multilogistic selected eliminating no useful data fields

modelfit_sel=stepAIC(glm(formula=MISSING ~.-FILARE-PIANTA-X-Y-candidate, family=binomial(), data=na.omit(mat_bulichella_sp@data)))
summary(modelfit_sel)

#############################################################################################################################
# model multilogistic NO selection but choices

formula_classifier_A1 ="MISSING ~ Area + Roughness + Line_res"
formula_classifier_A2 ="MISSING ~ Area + Roughness + Underperc"
formula_classifier_A2 ="MISSING ~ Area + Roughness + L_moran_50"

modelfit_A1 <- glm(formula=formula_classifier_A1 , family=binomial(), data=na.omit(mat_bulichella_sp@data))
summary(modelfit_A1)
modelfit_A2 <- glm(formula=formula_classifier_A2 , family=binomial(), data=na.omit(mat_bulichella_sp@data))
summary(modelfit_A2)
modelfit_A3 <- glm(formula=formula_classifier_A2 , family=binomial(), data=na.omit(mat_bulichella_sp@data))
summary(modelfit_A3)

sjt.glm(modelfit_A1 , modelfit_A2 , modelfit_A3 ,file="results/table_glm_compare.html")

sjt.glm(modelfit_A1,file="results/table_glm_A1.html")
sjt.glm(modelfit_A2,file="results/table_glm_A2.html")
sjt.glm(modelfit_A3,file="results/table_glm_A3.html")

##########################################################################################################

observed=mat_bulichella_sp@data$MISSING
prob_A1 <- predict(modelfit_A1, newdata=mat_bulichella_sp@data, type='response')
prob_A2 <- predict(modelfit_A2, newdata=mat_bulichella_sp@data, type='response')
prob_A3 <- predict(modelfit_A3, newdata=mat_bulichella_sp@data, type='response')

roc_data <- data.frame(Obs = observed, Model_A1 = prob_A1, Model_A3 = prob_A3, stringsAsFactors = FALSE)

longtest <- melt_roc(roc_data,"Obs", c("Model_A1","Model_A3"))

names(longtest)[3]="Classifier"
a=ggplot(longtest, aes(d = D, m = M, color = Classifier)) + geom_roc() + style_roc() 


a+annotate("text", x = .75, y = .25, 
           label = paste(" AUC A1=", round(calc_auc(a)$AUC[1], 2),"\n",
                         "AUC A3=", round(calc_auc(a)$AUC[2], 2))) 


ggsave("results/modelfit_ROC_glm.png",dpi=300)


res=list()
tresh_miss=c(0.05,0.1,0.15,0.2,0.8,0.9)
for ( i  in 1:length(tresh_miss)){
  prediction=ifelse(prob_A1>tresh_miss[i],1,0)
  res[[i]]=as.data.frame.array(table(observed,prediction))
}
res_df=do.call("rbind",lapply(res,function(x) data.frame(hit_plant=x[1,1],nohit_plants=x[1,2],nohit_miss=x[2,1],hit_miss=x[2,2])))
res_df$tresh=tresh_miss

sjt.df(res_df, describe = FALSE,show.rownames = FALSE, file = "results/model_results_glm_A1.html")

res=list()
tresh_miss=c(0.05,0.1,0.15,0.2,0.8,0.9)
for ( i  in 1:length(tresh_miss)){
  prediction=ifelse(prob_A2>tresh_miss[i],1,0)
  res[[i]]=as.data.frame.array(table(observed,prediction))
}
res_df=do.call("rbind",lapply(res,function(x) data.frame(hit_plant=x[1,1],nohit_plants=x[1,2],nohit_miss=x[2,1],hit_miss=x[2,2])))
res_df$tresh=tresh_miss

sjt.df(res_df, describe = FALSE,show.rownames = FALSE, file = "results/model_results_glm_A2.html")

res=list()
tresh_miss=c(0.05,0.1,0.15,0.2,0.8,0.9)
for ( i  in 1:length(tresh_miss)){
  prediction=ifelse(prob_A3>tresh_miss[i],1,0)
  res[[i]]=as.data.frame.array(table(observed,prediction))
}
res_df=do.call("rbind",lapply(res,function(x) data.frame(hit_plant=x[1,1],nohit_plants=x[1,2],nohit_miss=x[2,1],hit_miss=x[2,2])))
res_df$tresh=tresh_miss

sjt.df(res_df, describe = FALSE,show.rownames = FALSE, file = "results/model_results_glm_A3.html")


#############################################################################################################################
# decision trees models

treemodel_full_A1 <- rpart(formula(modelfit_A1), data=mat_bulichella_sp@data)

png("results/decision_tree_model_A1.png")
rpart.plot(treemodel_full_A1)
dev.off()

treemodel_full_A2 <- rpart(formula(modelfit_A2), data=mat_bulichella_sp@data)

png("results/decision_tree_model_A1.png")
rpart.plot(treemodel_full_A2)
dev.off()

treemodel_full_A3 <- rpart(formula(modelfit_A3), data=mat_bulichella_sp@data)

png("results/decision_tree_model_A3.png")
rpart.plot(treemodel_full_A3)
dev.off()

###############################################################################################################

prob_A1 <- predict(treemodel_full_A1, newdata=mat_bulichella_sp@data)
prob_A2 <- predict(treemodel_full_A2, newdata=mat_bulichella_sp@data)
prob_A3 <- predict(treemodel_full_A3, newdata=mat_bulichella_sp@data)

roc_data <- data.frame(Obs = observed, DTree_A1 = prob_A1, DTree_A3 = prob_A3, stringsAsFactors = FALSE)

longtest <- melt_roc(roc_data,"Obs", c("DTree_A1","DTree_A3"))

names(longtest)[3]="Classifier"
a=ggplot(longtest, aes(d = D, m = M, color = Classifier)) + geom_roc() + style_roc() 


a+annotate("text", x = .75, y = .25, 
           label = paste(" AUC A1=", round(calc_auc(a)$AUC[1], 2),"\n",
                         "AUC A3=", round(calc_auc(a)$AUC[2], 2))) 


ggsave("results/modelfit_ROC_tree.png",dpi=300)


#########################################################################################################################

res=list()
tresh_miss=c(0.05,0.1,0.15,0.2,0.8)
for ( i  in 1:length(tresh_miss)){
  prediction=ifelse(prob_A1>tresh_miss[i],1,0)
  res[[i]]=as.data.frame.array(table(observed,prediction))
}
res_df=do.call("rbind",lapply(res,function(x) data.frame(hit_plant=x[1,1],nohit_plants=x[1,2],nohit_miss=x[2,1],hit_miss=x[2,2])))
res_df$tresh=tresh_miss

sjt.df(res_df, describe = FALSE,show.rownames = FALSE, file = "results/model_results_tree_A1.html")

res=list()
tresh_miss=c(0.05,0.1,0.15,0.2,0.8)
for ( i  in 1:length(tresh_miss)){
  prediction=ifelse(prob_A2>tresh_miss[i],1,0)
  res[[i]]=as.data.frame.array(table(observed,prediction))
}
res_df=do.call("rbind",lapply(res,function(x) data.frame(hit_plant=x[1,1],nohit_plants=x[1,2],nohit_miss=x[2,1],hit_miss=x[2,2])))
res_df$tresh=tresh_miss

sjt.df(res_df, describe = FALSE,show.rownames = FALSE, file = "results/model_results_tree_A2.html")

res=list()
tresh_miss=c(0.05,0.1,0.15,0.2,0.8)
for ( i  in 1:length(tresh_miss)){
  prediction=ifelse(prob_A3>tresh_miss[i],1,0)
  res[[i]]=as.data.frame.array(table(observed,prediction))
}
res_df=do.call("rbind",lapply(res,function(x) data.frame(hit_plant=x[1,1],nohit_plants=x[1,2],nohit_miss=x[2,1],hit_miss=x[2,2])))
res_df$tresh=tresh_miss

sjt.df(res_df, describe = FALSE,show.rownames = FALSE, file = "results/model_results_tree_A3.html")





# A. Baddeley, E. Rubak and R.Turner. Spatial Point Patterns: Methodology and Applications with R. Chapman and Hall/CRC Press, 2015.