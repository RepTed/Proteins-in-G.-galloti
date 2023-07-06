#Environmental factors influence biomarkers in the lizard Gallotia galloti
#------------------------Microclimate------------------------

##=======Temperature variables=======

#A seperate script was written to better estimate temperature on the microclimate scale

install.packages("devtools")
library(devtools)
install_github('mrke/NicheMapR')
devtools::install_github('ilyamaclean/microclima')
library(raster)
library(NicheMapR)
library(microclima)


#Site:BV
#Input coordinates here for site location. Resolution is in meters. 
loc <- microclima::get_dem(lat = 28.36816, long = -16.86408, resolution = 30)

#Run the "habitats" line to look at the available habitat options
habitats
#Using habitat 15 (Cropland/Natural vegetation mosaic) for BV, TA, ST, CR and LL.
#Using habitat 16 (Barren or sparsely vegetated) for TSE and TRE.

#Input into following function. Date is 17th July 2018, at 3cm above ground.
temps.bv <- runauto(loc, "17/07/2018", "17/07/2018", hgt = 0.03,
                    l = NA, x = NA,
                    habitat = "Cropland/Natural vegetation mosaic",
                    plot.progress = FALSE)

?runauto
#Visualise the mean temperatures across the loc raster by:
plot(temps.bv$tmean)

#Extract the relevant information:
#Convert mean temperature list to data frame.
temp.df <- as.data.frame(temps.bv[["tmean"]]@data@values)
#Remove any NA's - mostly sea areas. 
temp.df <- na.omit(temp.df)
#Calculate mean for the location.
bv.mean <- round(mean(temp.df$`temps.bv[["tmean"]]@data@values`), digits = 2)
#Extract max temperature
bv.max <- round(temps.bv[["tmax"]]@data@max, digits = 2)
#Extract min temperature
bv.min <- round(temps.bv[["tmin"]]@data@min, digits = 2)
#Add site name
bv <- "BV"
#Combine into data frame
bv.df <- as.data.frame(cbind(bv, bv.mean, bv.min, bv.max))
#Change column names to easily combine later
colnames(bv.df) <- c("Site", "Mean", "Min", "Max")

#Repeat for each locality by changing date and coordinates


#Final dataset
#Combine all data frames
final.data <- rbind(bv.df, ta.df, st.df, cr.df, ll.df, tse.df, tre.df)
#Save to a csv file - need to input new filepath after file = 
library(readr)
write_csv(file = "Microclimate_sites.csv", x = final.data)


##=======script for generating microclimate variables=======

#library(devtools)
#install_github('mrke/NicheMapR')
#devtools::install_github('ilyamaclean/microclima')
library(raster)
library(NicheMapR)
library(microclima)

#Site:BV
#Run model
micro.bv <- micro_ncep(loc = c(28.36816, -16.86408), "17/07/2018", "17/07/2018",
                       dem.res = 30, Usrhyt = 0.03)

##metout info is for minimum specified shade##
#Change relevant output to data frame
metout.bv <- as.data.frame(micro.bv$metout)
#Extract air temp at 3cm
bv.taloc <- mean(metout.bv$TALOC)
#Extract air temp at 2m
bv.taref <- mean(metout.bv$TAREF)
#Extract relative humidity at 3cm
bv.rhloc <- mean(metout.bv$RHLOC)
#Extract relative humidity at 2m
bv.rh <- mean(metout.bv$RH)
#Extract wind speed at 3cm
bv.vloc <- mean(metout.bv$VLOC)
#Extract wind speed at 2m
bv.vref <- mean(metout.bv$VREF)
#Extract soil surface wetness
bv.pctwet <- mean(metout.bv$PCTWET)
#Extract zenith angle of sun
bv.zen <- mean(metout.bv$ZEN)
#Extract solar radiation
bv.radi <- mean(metout.bv$SOLR)
#Extract sky radiant temperature
bv.tskyc <- mean(metout.bv$TSKYC)


#add site name too
site <- "BV"

#combine values
metout.data <- as.data.frame(cbind(site, bv.taloc, bv.taref, bv.rhloc, bv.rh, bv.vloc, bv.vref,
                                   bv.pctwet, bv.zen, bv.radi, bv.tskyc))

#rename columns to easily stick new data frames together
colnames(metout.data) <- c("Site", "TALOC", "TAREF", "RHLOC", "RH", "VLOC", 
                           "VREF", "PCTWET", "ZEN", "SOLR", "TSKYC")

#combine data sets
bv.data <- metout.data

#Repeat for each locality by changing date and coordinates

#Final dataset
#Combine all data frames
final.data <- rbind(bv.data, ta.data, st.data, cr.data, ll.data,
                    tse.data, tre.data)
#Save to a csv file - need to input new filepath after file = 
library(readr)
write_csv(file = "Final_data_v2.csv", x = final.data)

###------------------------model selection------------------------

#Generalized linear model selection using GLMulti

library(glmulti)
library(dplyr)
library(flextable)

mydata <- read.csv("data.csv")
attach(mydata)

#Check R knows these are factor variables not integers
mydata$Env <- as.factor(mydata$Env)
mydata$Morphotype <- as.factor(mydata$Morphotype)
mydata$Sex <- as.factor(mydata$Sex)
mydata$Population <- as.factor(mydata$Population)

#A model using all the variables included is below for reference
allmod <- glm(GRP94p_TPF ~Morphotype + Env + Sex + Elevation + Daily_ave_.temp + Population + Aspect_val + Slope_val +
                 Proc_time + RH + PCTWET +SOLR + TSKYC, family = Gamma(link ="log"), data = mydata)

#g = genetic algorithm
#h = exhaustive run through all models
#d = number of models

#run exhaustive combination of all models WITHOUT interactions
h_model <- glmulti(GRP94p_TPF ~ Morphotype + Env + Sex + Elevation + Daily_ave_.temp + Population + Aspect_val + Slope_val +
                     Proc_time + RH + PCTWET +SOLR + TSKYC,
                   data = mydata,
                   crit = aicc, 
                   level = 1, 
                   method = 'h',
                   family = Gamma(link ="log"), 
                   fitfunction = glm, 
                   confsetsize = 100)

#compare the weights and aicc of top X amount models, that are within deltaAICc of 2
weightable(h_model)[1:5,] %>% 
  regulartable() %>% 
  autofit()

#Plot the most frequently occurring model-averaged importance of terms
plot(h_model, type = "s")
#The 80% cut off is somewhat arbitrary, I included the top four, in this case they also appear in the top models

#what variables are included in all models with deltaAICc of 2?
#Slope_val
#SOLR
#Daily_ave_.temp
#TSKYC
#RH
#Elevation

#Re-run exhaustive model selection including above terms, WITH interactions (level = 2)
h_model <- glmulti(GRP94p_TPF ~ Elevation + Daily_ave_.temp + Slope_val +
                     RH + SOLR + TSKYC,
                   data = mydata,
                   crit = aicc, 
                   level = 2, 
                   method = 'h',
                   family = Gamma(link ="log"), 
                   fitfunction = glm, 
                   confsetsize = 100)

#Repeat above validation techniques, looking at all the models with equal support (deltaAICc = 2), and top four model-averaged important terms
#Top models with equal support
weightable(h_model)[1:5,] %>% 
  regulartable() %>% 
  autofit()
#model-averaged importance of terms
plot(h_model, type = "s")
#Compare which terms are included from the table and the plot, and use these to inform which predictors to include in the final model, consdering the final AICc score
#in this case the final model looks like this:
GRP94mod <- glm(GRP94p_TPF ~ Slope_val:Elevation + TSKYC:RH, family = Gamma(link ="log"), data = mydata)

#Repeat this for all biomarkers
#family=Gamma(link="log") was used for all except total 3-NT, which was family = gaussian

##other useful tool for checking models

library(performance)

#retrieve information on model
model_performance(GRP94mod)
#check and visualise all assumptions of model
check_model(GRP94mod)
library(effects)
#visualising model
plot(allEffects(GRP94mod), lines = list(multiline = T), confint = list(style = "bars"))

###------------------------Collinearity plots------------------------

library(corrplot)
#You need a dataframe with just the variables you are interesting in testing for collinearity
mydata <- read.csv("Predictors_only.csv")

M<-cor(mydata)
corrplot(M, addCoef.col = "black",number.cex=0.4)

#row and columns can then be labelled with the full variable names as appropriate
colnames(M) <- c("1",
                 "2")
rownames(M) <- c("1",
                 "2")
#other useful tool for checking collinearity, build a model with all variables in and then;
library(performance)
check_collinearity(all.variables.mod)
plot(check_collinearity(all.variables.mod))

###------------------------Plot Interactions and post-hoc emmeans------------------------
library(interactions)
library(scales)
library(ggplot2)
#Plotting the interaction plots, with partial residues between the significant interaction terms from each of the final models
#make glms
carbtot <- glm(carbtot ~ Sex + TSKYC + SOLR:Elevation, family = Gamma(link ="log"), data = mydata)
grp94 <- glm(GRP94p_TPF ~ Slope_val:Elevation+TSKYC:RH, family = Gamma(link ="log"), data = mydata)
#total carbonylation interaction
carbtotint<-interact_plot(carbtot, pred = Elevation, modx = SOLR, plot.points = TRUE, partial.residuals = T, 
                          interval = T, int.width = 0.5, point.size = 2, 
                          y.label = "Relative total carbonylation (A.U)", x.label = "Elevation (m)") + scale_y_continuous(trans = "log10",
                                                                                                                          labels = trans_format("log10", math_format(10^.x))) + theme_classic(base_size = 16) + expand_limits(y=c(10^4, 10^5))
#grp94 interaction 1
grpint1<- interact_plot(grp94, pred = Slope_val, modx = Elevation, plot.points = TRUE, partial.residuals = T, 
                        interval = T, int.width = 0.8, point.size = 2, x.label = "Slope value (arc °)", 
                        y.label = "Relative GRP94 expression (A.U)") + scale_y_continuous(trans = "log10",
                                                                                          labels = trans_format("log10", math_format(10^.x))) + theme_classic(base_size = 16) + expand_limits(y=c(10^3, 10^5))

#grp94 interaction 2
grpint2<-interact_plot(grp94, pred = TSKYC, modx = RH, plot.points = TRUE, partial.residuals = T, 
                       interval = T, int.width = 0.8, point.size = 2, x.label = "Sky radiant temperature (°C)", 
                       y.label = "Relative GRP94 expression (A.U)") + scale_y_continuous(trans = "log10",
                                                                                         labels = trans_format("log10", math_format(10^.x))) + theme_classic(base_size = 16) + expand_limits(y=c(10^3, 10^5))
#post-hoc analysis for categorical variable
library(emmeans)
emmdf <- read.csv("sex_emm.csv")
sexemm <- ggplot(emmdf, aes(x = Sex, y = emmean)) + theme_classic(base_size = 16) +
  geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width = 0.1, colour = "cornflowerblue", linewidth = 2) +
  geom_point(size = 5, colour = "lightblue") + ylab("Estimated marginal means \n (Relative total carbonylation)")
###------------------------Modelling Te------------------------
#Modelling the operative temperature
library(TrenchR)
#Detailed instructions found here: https://cran.r-project.org/web/packages/TrenchR/vignettes/TeTutorial.html
#Section: Example operative environmental temperature calculations
#Using environmental inputs from microclima model above
#Using model from Fei

