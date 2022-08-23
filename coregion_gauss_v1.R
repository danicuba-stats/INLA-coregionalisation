################################################################################
##### Co-Regionalisation model for Heavy Metal pairs in the Glasgow region #####
##### Date: 28th of July, 2022
##### Author: Daniela Cuba
################################################################################

# Load some libs
libs<-c("sp","maptools","dplyr","geoR","gstat",
        "INLA","rgdal","broom","gridExtra","raster","gtools",
        "ggplot2")
lapply(libs,require,character.only=T)

# Load workspace
load("C:/Users/student/OneDrive - University of Glasgow/Documents/PhD/GBASE_INLA+excursion/coregionalisation/coregion_gauss_v1_env.RData")

# If not, run everything here 
#### Data ----
# Read G-BASE data
data<-read.csv("C:/Users/student/OneDrive - University of Glasgow/Documents/PhD/GBASE_data_source/gbase_glasgow_covs.csv")
# Read prediction covariates
data_pred<-read.csv("C:/Users/student/OneDrive - University of Glasgow/Documents/PhD/GBASE_data_source/gbase_prediction_covs.csv")

# Subset desireable variables
elems<-c("As","Cr","Cu","Pb","Ni","Zn")
data<-data[,c("X_COORD","Y_COORD",
              "As_XRF","Cr_XRF","Cu_XRF","Pb_XRF","Ni_XRF","Zn_XRF",
              "elevation","slope","aspect","plan.curve","profile.curve",
              "twi","mrvbf","mrrtf","population","landuse")]
data$Pb_XRF<-as.numeric(gsub("[^0-9\\.]", "", data$Pb_XRF))
data$Cr_XRF<-as.numeric(gsub("[^0-9\\.]", "", data$Cr_XRF))
data$As_XRF<-as.numeric(gsub("[^0-9\\.]", "", data$As_XRF))
data$Zn_XRF<-as.numeric(gsub("[^0-9\\.]", "", data$Zn_XRF))
data$Cu_XRF<-as.numeric(gsub("[^0-9\\.]", "", data$Cu_XRF))
data$Ni_XRF<-as.numeric(gsub("[^0-9\\.]", "", data$Ni_XRF))

data<-data[complete.cases(data),]

names(data)[1:8]<-c("X","Y",elems)

data.sp<-SpatialPoints(coords=data[,c("X","Y")])

data$logAs<-log(data$As)
data$logCr<-log(data$Cr)
data$logCu<-log(data$Cu)
data$logPb<-log(data$Pb)
data$logNi<-log(data$Ni)
data$logZn<-log(data$Zn)

# Prepare element pairs
elem_pairs<-t(apply(combinations(6,2,1:6),1,function(x) elems[x]))
pair_names<-vector("character",15)
for(i in 1:nrow(elem_pairs)){
  pair_names[i]<-paste(elem_pairs[i,1],elem_pairs[i,2],sep="-")
}

#### Spatial data ----
# Spatial information
counties <- readOGR("C:/Users/student/OneDrive - University of Glasgow/Documents/Other Data/UK_county_boundaries/CountiesOS.shp")
data.sp <- SpatialPointsDataFrame(coords = data[,c("X","Y")],
                                  data=data.frame(logCr=data$logCr,
                                                  Cr=data$Cr),
                                  proj4string = CRS(proj4string(counties)))
gbase_counties <- over(data.sp, counties)
gbase_counties <- gbase_counties$NAME %>% table
gbase.pols <- counties[which(counties$NAME%in%names(gbase_counties)),]

gbase.counties.data <- tidy(gbase.pols) 
# names(gbase.counties.data)[1:2] <- c("Eastings","Northings")

#### Create Model ----
# Define mesh
mesh <- inla.mesh.2d(loc=data[,c("X","Y")], 
                     loc.domain=data_pred[,c("X","Y")],
                     max.edge = c(7500,30000),
                     offset = c(5000,20000))
# plot(mesh)

# Build SPDE
spde<-inla.spde2.pcmatern(mesh=mesh,
                          alpha=2,
                          prior.range=c(10000,0.05),
                          prior.sigma=c(2,0.05)) # Priors that worked well from the sensitivity analysis
# At the moment, priors are the same for all pairs
A.pred<- inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(data_pred[,c("X","Y")]))
A.obs<- inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(data[,c("X","Y")]))
s.index<-inla.spde.make.index(name="spatial.field",
                              n.spde=spde$n.spde)

# Hyperparameters
hyper_norm <- list(beta = list(prior = 'normal', param = c(0, 10)))
hyper_eps<- list(hyper=list(prec=list(prior = "pc.prec", 
                                      param=c(1,0.01))))

# Formula - I am not confident this is correct
form <- y ~ 0 + intercept1 + intercept2 + f(s1,model=spde) +
  f(s2, copy="s1", fixed=F, hyper=hyper_norm)

# Create inla stack for observed data
stack1<-inla.stack(data=list(y=cbind(as.vector(data$logCr),NA)),
                  A=list(A.obs),
                  effects = list(list(s1=s.index,intercept1=1)), tag="logCr.data")
stack2<-inla.stack(data=list(y=cbind(NA,as.vector(data$logAs))),
                  A=list(A.obs),
                  effects = list(list(s2=s.index,intercept2=1)), tag="logAs.data")
stack3<-inla.stack(data=list(y=cbind(rep(NA,nrow(data_pred)),
                                     rep(NA,nrow(data_pred)))),
                   A=list(A.pred),
                   effects = list(list(s2=s.index,intercept2=1)), tag="pred")
monster.stack<-inla.stack(stack1,stack2,stack3)


# Fit INLA model
m.logCr<-inla(form,
              data=inla.stack.data(monster.stack,
                                   spde=spde),
              family=c("gaussian","gaussian"),
              #control.family=list(hyper_norm, hyper_norm),
              control.predictor=list(A=inla.stack.A(monster.stack), compute=TRUE),
              #control.inla=list(theta=theta.ini,restart=T),
              control.compute=list(cpo=TRUE,dic=TRUE,config=TRUE))
summary(m.logCr)

id.prd <- inla.stack.index(monster.stack, "pred")$data # obtain the index of the points to predict

# Extract marginals (logCr and logAs)
pred.vars<- do.call(rbind,lapply(m.logCr$marginals.fitted.values,function(x) apply(x,2,mean)))
new.pred.grid_gaus<-as.data.frame(data_pred)
new.pred.grid_gaus$pred<- m.logCr$summary.fitted.values$mean[id.prd]
new.pred.grid_gaus$logCr<-pred.vars[id.prd,1]
new.pred.grid_gaus$logAs<-pred.vars[id.prd,2]

## Plot Cr
ggplot() +  
  geom_point(data = new.pred.grid_gaus,
             aes(x = X, y = Y,col=logCr),size=2.1,shape=15) +
  geom_polygon(data=gbase.counties.data,
               aes(x = long, y = lat, group = group),
               col = "black", fill = NA, size=1.2) +
  coord_cartesian(xlim = range(data$X),
                  ylim = range(data$Y)) +
  # scale_color_viridis_c(limits=range(c(new.pred.grid$pred,new.pred.grid_gaus$pred)))+
  scale_color_viridis_c()+
  labs(x = "X",y = "Y", legend = "logCr") +#+
  ggtitle("Gaussian Model - Marginal:logCr")#+

## Plot As
ggplot() +  
  geom_point(data = new.pred.grid_gaus,
             aes(x = X, y = Y,col=logAs),size=2.1,shape=15) +
  geom_polygon(data=gbase.counties.data,
               aes(x = long, y = lat, group = group),
               col = "black", fill = NA, size=1.2) +
  coord_cartesian(xlim = range(data$X),
                  ylim = range(data$Y)) +
  # scale_color_viridis_c(limits=range(c(new.pred.grid$pred,new.pred.grid_gaus$pred)))+
  scale_color_viridis_c()+
  labs(x = "X",y = "Y", legend = "LogAs") +#+
  ggtitle("Gaussian Model - Marginal:LogAs")#+


##### PROBLEM: 
new.pred.grid_gaus$logCr %>% range
# this range is 1.883903 to 2.758874
# but the actual range of the data is 
data$logCr %>% range
# 3.16 to 8.58. 
# Either the model is incorrectly specified, or I am extracting the incorrect information. 
# This is also the case with the second variable in the model, LogAs
new.pred.grid_gaus$logAs %>% range
# this range is between 1.46 and 4.34
# but the actual range of the variable is
data$logAs %>% range
# -0.35667 to 6.7452
#########################################################################################


# save.image("C:/Users/student/OneDrive - University of Glasgow/Documents/PhD/GBASE_INLA+excursion/coregionalisation/coregion_gauss_v1_env.RData")
# TEst - curious but it doesn't work
# library(excursions)
# res.qc_gaus<-excursions.inla(m.logCr,stack=monster.stack,
#                              tag="pred", alpha=0.99,u=log(100),method='QC',type=">",max.threads=0)
##### Doesn't work - everything is 0 or NA
