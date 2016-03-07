# Individual recognition and missing plant detection in remotely sensed vineyard aerial maps. 

**Authors:** 

*Jacopo Primicerio [1,2] , Giovanni Caruso [3], Lorenzo Comba [2], Alfonso Crisci [1], Paolo Gay [2,4], Lorenzo Genesio [1], Francesco Primo Vaccari [1]*

**Author's affiliation:** 

**[1] CNR – IBIMET**: Istituto di Biometeorologia  Via G.Caproni, 8, 50145 – Italy, j.primicerio@ibimet.cnr.it, a.crisci@ibimet.cnr.it, l.genesio@ibimet.cnr.it, f.vaccari@ibimet.cnr.it 

**[2] D.A.F.E.**: Università degli Studi di Pisa, 80 Via del Borghetto, 56124 Pisa – Italy giovanni.caruso@for.unipi.it

**[3] DI.S.A.F.A.**:  Università degli Studi di Torino, 44 Via Leonardo da Vinci, 10095 Grugliasco (TO) – Italy lorenzo.comba@unito.it, paolo.gay@unito.it

**[4] CNR-I.E.I.I.T**: Istituto di IEIIT, Corso Duca degli Abruzzi  24 , 10129 Torino – Italy 

**Corresponding author:** 

*Jacopo Primicerio* Email: j.primicerio@ibimet.cnr.it Tel: +39 055 3033711


## Material and methods

### Missing plants detection

Previous tasks of the work has given a set of spatial features, with own specific attributes, that are located along the lines of the vineyard. Teorethically each one represent the signature on the UAV image of single plant. Feature identification was built by geometrical costruction. A former data matrix is available and give the plant centroid coordinates ( X and Y) , the  correpondent area and perimeter of feature and lines belonging. Starting from these data the core aim of the work is to detect which are the features that to be considered as missing plants and if possible to design a spatial pattern of missingness on the field. Strictly speaking the general hypothesys is to verify if features are real plants or not. The simplest way to reach this target is to build a statistical classifier  based on a multiple logistic regression scheme. Logistic regression works with independent variables, named predictors, and give a probabilistic prediction concerning  a dependent variable that could be assimilated a binary categorical state. In our cases become  a boolean one : *missing plant* (yes /not). Following the way a  generalized linear model (GLM) procedure  was implemented in R Stat environment ( R Core Team 2015). A comprehensive training set was built considering all data available. A boolean predictand variable, named as **MISSING**, was built mapping which are features associated with plant missingness. Several predictors variables are considered. The first one used, the main in order of importance, is the **Area** of the plant's feature. Seems very reasonable  arguing that low area values corresponds to the voids left in the line by plant's disappearance. This parameter has deeply investigated and a plot of quantile distribution was done in figure X. It is hard to define a fixed treshshold in the area parameter to identificate missing plants because plant's vigour is not spatially uniform in vineyards. Neighbours plant in function to their vigour produce a relevant bias on area parameter reducing the skill of this predictor.The second one consists in the **Perimetro**  and consist in the length of polygonal boundary of features. This predictor obviously is connected to area and in model but could show collinearity. The third is a boolean factor predictor, named **Underperc** and inform if area's feature is below to 10th percentile area treshshold. This predictor is a information enhancer of plant missingness taking into account the comprehensive distribution of area parameter. Last tree predictors, named respectively as **L_moran**,**dist_bound** and **Roughness**, are different to previous ones and are calculated thanks to specific geoprocessing procedures applied on data matrix converted in GIS point object.Concerning the first two predictors the R packages used are: *raster* (Hijmans,2015)and *rgeos* (Bivand and Colin Rundel,2015).For the last one other GIS procedures was carried out. **L_moran** predictor is the value of local Moran's I index calculated for each feature and give a measure of spatial association concerning area paraemeter. It is known that plant missingness could be associated to Moran's I minimum (Filzmoser et al., 2014) because missing plants showing patterns with spatial density maximum where low plant area values show continuity.**dist_bound** is another geometrical predictor taking into account that measures the distance of plants from the field boundary. It is used to consider the eventual bound-effect of vineyard on plant missingness. Last index calculated and used as predictor is **Roughness** that is a measure of “compactness” of feature following  referenced methodologies (Tang and Tian, 2008). Substantially it's a parameter shape dependent  and gives information in regard to the morphology of plant soil covering. It could be useful to measure influence of the bias due to neighbour plants especially when plant missingness occurs. In order to individuate the best model, hence the more perfomant classifier, a model selection was carried out by using area under curve criteria (AUC) of receiver operating characteristic (ROC)  of outcomes face to field training data. Starting from a full model, where all predictors, the  useless one were discarded. Afterwise three different propability level was analized ( p =0.95,0.9,0.85) a categorical model response carried out on  best model previously selected. At this regard Brier Score Index of multi logistic classifier was done by using  R *verification* (NCAR, 2014) and *ROCR* package (Sing et al,2005).

## Results

Full multiple logistic model statistics  are summarized in the following table. Odd ratio is the common way to  describes the association between the presence/absence of dependent with the predictor/s. **Area** and **L_moran** (Local Moran's I) confirms, as we expect, their strong skill for plant missing detection showing  great significance. A second group of predictors ( **Underperc** and **Roughness**) show a weak significance indicate that are useful only as auxilllary variables. Other ones ( **Perimetro** and **dist_bound**) doesn't show any influence on missingness detection. 

**Table 1** *Summary of complete multi logistic model.*

<table style="border-collapse:collapse; border:none;border-bottom:double;"><tr>
<td style="padding:0.2cm; border-top:double;">&nbsp;</td>
<td style="padding-left:0.5em; padding-right:0.5em; border-top:double;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; border-bottom:1px solid; border-top:double;" colspan="3">Model complete</td>
</tr>
<tr>
<td style="padding:0.2cm; font-style:italic;">&nbsp;</td><td style="padding-left:0.5em; padding-right:0.5em; font-style:italic;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; font-style:italic; ">Odd Ratios</td><td style="padding:0.2cm; text-align:center; font-style:italic; ">CI</td><td style="padding:0.2cm; text-align:center; font-style:italic; ">p</td> 
</tr>
<tr>
<td style="padding:0.2cm; border-top:2px solid; text-align:left;">Intercept</td><td style="padding-left:0.5em; padding-right:0.5em; border-top:2px solid; ">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; border-top:2px solid; ">150.63</td><td style="padding:0.2cm; text-align:center; border-top:2px solid; ">0.10&nbsp;&ndash;&nbsp;27168.29</td><td style="padding:0.2cm; text-align:center; border-top:2px solid; ">.098 . </td>
</tr>
<tr>
<td style="padding:0.2cm; text-align:left;">Area</td><td style="padding-left:0.5em; padding-right:0.5em;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; ">0.99</td><td style="padding:0.2cm; text-align:center; ">0.98&nbsp;&ndash;&nbsp;0.99</td><td style="padding:0.2cm; text-align:center; "><b>.001 ***</b></td>
</tr>
<tr>
<td style="padding:0.2cm; text-align:left;">Perimetro</td><td style="padding-left:0.5em; padding-right:0.5em;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; ">1.01</td><td style="padding:0.2cm; text-align:center; ">0.96&nbsp;&ndash;&nbsp;1.10</td><td style="padding:0.2cm; text-align:center; ">.698</td>
</tr>
<tr>
<td style="padding:0.2cm; text-align:left;">L_moran</td><td style="padding-left:0.5em; padding-right:0.5em;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; ">1.16</td><td style="padding:0.2cm; text-align:center; ">1.00&nbsp;&ndash;&nbsp;1.25</td><td style="padding:0.2cm; text-align:center; "><b>.002 **</b></td>
</tr>
<tr>
<td style="padding:0.2cm; text-align:left;">dist_bound</td><td style="padding-left:0.5em; padding-right:0.5em;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; ">0.00</td><td style="padding:0.2cm; text-align:center; ">0.00&nbsp;&ndash;&nbsp;Inf</td><td style="padding:0.2cm; text-align:center; ">.291</td>
</tr>
<tr>
<td style="padding:0.2cm; text-align:left;">Underperc</td><td style="padding-left:0.5em; padding-right:0.5em;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; ">0.51</td><td style="padding:0.2cm; text-align:center; ">0.24&nbsp;&ndash;&nbsp;1.09</td><td style="padding:0.2cm; text-align:center; ">.085 .</td>
</tr>
<tr>
<td style="padding:0.2cm; text-align:left;">Roughness</td><td style="padding-left:0.5em; padding-right:0.5em;">&nbsp;</td>
<td style="padding:0.2cm; text-align:center; ">0.01</td><td style="padding:0.2cm; text-align:center; ">0.00&nbsp;&ndash;&nbsp;4.64</td><td style="padding:0.2cm; text-align:center; ">.081 .</td>
</tr>
<tr>
<td style="padding:0.2cm; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;">N° Observations</td><td style="padding-left:0.5em; padding-right:0.5em; border-top:1px solid;">&nbsp;</td> <td style="padding:0.2cm; padding-top:0.1cm; padding-bottom:0.1cm; text-align:center; border-top:1px solid;" colspan="3">2242</td>
</tr>
<tr>
<td style="padding:0.2cm; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">AIC value </td><td style="padding-left:0.5em; padding-right:0.5em;">&nbsp;</td><td style="padding:0.2cm; text-align:center; padding-top:0.1cm; padding-bottom:0.1cm;" colspan="3">728.904</td>
</tr>
</table>



##References

[@R01]: R  Core Team (2015). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.  <URL: http://www.R-project.org/>.

[@R02]:Robert J. Hijmans (2015). raster: Geographic Data Analysis and Modeling. R package version 2.4-18.  <URL:http://CRAN.R-project.org/package=raster>.

[@R03]:Roger Bivand and Colin Rundel (2015). rgeos: Interface to Geometry Engine - Open Source (GEOS). R package version 0.3-12.  <URL:http://CRAN.R-project.org/package=rgeos>.

[@R04]:NCAR - Research Applications Laboratory (2014). verification: Weather Forecast Verification Utilities.. R package version 1.41. http://CRAN.R-project.org/package=verification
  
[@R05]: Sing T, Sander O, Beerenwinkel N and Lengauer T (2005). “ROCR: visualizing classifier performance in R.” _Bioinformatics_, *21*(20), pp. 7881. <URL:http://rocr.bioinf.mpi-sb.mpg.de>.  
  
**Contacts**:

IBIMET CNR Jacopo Primicerio Emails : copojax@gmail.com , j.primicerio@ibimet.cnr.it
