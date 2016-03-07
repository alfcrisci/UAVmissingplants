#############################################################
# Setup work directory

setwd("/home/alf/Scrivania/documents/lav_primicerio/final/code")

#############################################################
# Load dependencies and functions

source("load_libraries.r")
source("auxillary_bulichella.r")

################################################################################################
# load geo data
setwd("/home/alf/Scrivania/documents/lav_primicerio/final/")

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


# Dist_bound = Calcolo della forma convessa e della distanza

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

# Quantile UnderPerc Area quantile sotto il 90 percentile variabile booleana

q_area=quantile(values(r_bulichella[["Area"]]),probs = seq(0, 1, 0.1),na.rm=T)
mat_bulichella_sp$Underperc=extract(r_bulichella[["Area"]]>q_area[10],mat_bulichella_sp)


####################################################################################

mat_bulichella_miss=mat_bulichella_sp[mat_bulichella_sp$MISSING>0,]

##################################################################################
# Create spatial object of missing density


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
# Local moran calculation. create matrix of neighbours


mat_bulichella_5=nb2listw(knn2nb(knearneigh(mat_bulichella_sp,k=5))) # five plants
mat_bulichella_100=nb2listw(knn2nb(knearneigh(mat_bulichella_sp,k=100))) # 100 plants
mat_bulichella_300=nb2listw(knn2nb(knearneigh(mat_bulichella_sp,k=300))) # 300 plants


mat_bulichella_sp@data$L_moran_5=localmoran(mat_bulichella_sp$Area, mat_bulichella_5)[,1]
mat_bulichella_sp@data$L_moran_100=localmoran(mat_bulichella_sp$Area, mat_bulichella_100)[,1]




###############################################################################################
# line vigour model model_primicerio function.
mat_bulichella_ls=split(mat_bulichella_sp@data,mat_bulichella_sp@data$FILARE)

res=list()
for (i in seq_along(mat_bulichella_ls)) {
  res[[i]]=model_primicerio(mat_bulichella_ls[[i]],
                            saveplot=TRUE,
                            titlefig=paste0("Modeling Plant Missing FILARE_",as.character(i)),
                            namefig=paste0("results/Modeling_Plant_Missing_FILARE ",as.character(i),".png"),
                            treshshold=150)
}


ls_model_residuals=lapply(res,function(x) x$model_residuals)

mat_bulichella_sp@data$Line_res=unlist(ls_model_residuals)

ls_canditate=lapply(res,function(x) x$vector)

mat_bulichella_sp@data$candidate=unlist(ls_canditate)

names(mat_bulichella_sp@data)


# [1] "FILARE"       "PIANTA"       "X"            "Y"            "Area"         "Perimetro"   
# [7] "Roughness"    "MISSING"      "LAT"          "LON"          "dist_bound"   "L_moran_30"  
# [13] "L_moran_50"   "L_moran_30_p" "L_moran_50_r" "Underperc"    "RatioAP"      "candidate"   


#############################################################################################################################
# model multilogistic


formula_classifier ="MISSING ~ Area + Roughness + Underperc + L_moran_50_r + Line_res"

modelfit_Asel <- glm(formula=formula_classifier , family=binomial(), data=na.omit(mat_bulichella_sp@data))
summary(modelfit_Asel)

##########################################################################################################
require(sjPlot)
lab <- c("Area", "Roughness", "Underperc", "L_moran_50_r","Line_res")
labdep <- c("Missingness")
sjt.glm(modelfit_Asel,
        labelDependentVariables=labdep,
        labelPredictors=lab,
        file="table_glm_predictors.html")


prob <- predict(modelfit_Asel, newdata=mat_bulichella_sp@data, type='response')
observed=mat_bulichella_sp@data$MISSING

res=list()
tresh_miss=c(0.05,0.1,0.15,0.2,0.8,0.9)

for ( i  in 1:length(tresh_miss)){
prediction_8=ifelse(prob>tresh_miss[i],1,0)
res[[i]]=as.data.frame.array(table(observed,prediction_8))
}

res_df=do.call("rbind",lapply(res,function(x) data.frame(hit_plant=x[1,1],nohit_plants=x[1,2],nohit_miss=x[2,1],hit_miss=x[2,2])))
res_df$tresh=tresh_miss
sjt.df(res_df, describe = FALSE, file = "model_results.html")

rocplot(modelfit_Asel)+ggtitle("Model GLM :Area + L_moran_R_50+ Roughness + Line_model_residuals ")
ggsave("results/modelfit_A.png",dpi=300)

saveRDS(modelfit_Asel,"modelfit_Asel.rds")
