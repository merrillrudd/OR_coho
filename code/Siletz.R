rm(list=ls())

####################
## load packages
####################

# devtools::install_github("merrillrudd/FishStatsUtils", ref="stream")
# devtools::install_github("james-thorson/VAST", ref='development')
# devtools::install_github("merrillrudd/VASTPlotUtils")

library(VAST)
library(TMB)
library(tidyverse)
library(VASTPlotUtils)

###################
## Directories
###################

main_dir <- "C:\\merrill\\OR_coho"
data_dir <- file.path(main_dir, "data")

sil_dir <- file.path(main_dir, "Siletz")
dir.create(sil_dir, showWarnings=FALSE)

fig_dir <- file.path(sil_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

####################
## Read in data
####################

network_all <- readRDS(file.path(data_dir, "network.rds"))
obs_all <- readRDS(file.path(data_dir, "observations.rds"))
hab_all <- readRDS(file.path(data_dir, "habitat.rds"))

## subset Siletz
network <- readRDS(file.path(data_dir, "network_siletz_v1.rds")) %>% rename(Lat = lat, Lon = long)
obs <- obs_all %>% filter(Population == "Siletz") %>% filter(dist_i > 0)
hab <- hab_all %>% filter(Population == "Siletz")

####################
## Network
####################
net_nodes <- network$child_s
Network_sz_LL <- network
Network_sz <- network %>% select("parent_s", "child_s", "dist_s")
plot_network(Network_sz_LL = Network_sz_LL, arrows=TRUE, root=TRUE)
plot_network(Network_sz_LL = Network_sz_LL, arrows=TRUE, root=TRUE, FilePath=fig_dir)

####################
## Observations
####################
obs_nodes <- obs$child_i
all(obs_nodes %in% net_nodes)

Data_count <- data.frame( "Catch_KG" = obs$Count, 
              "Year" = as.numeric(obs$Year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = obs$dist_i,
               "Lat" = obs$Lat, 
               "Lon" = obs$Lon, 
               "Pass" = 0,
               "Knot" = unlist(obs$child_i),
               "Category" = obs$Survey,
               "CategoryNum" = sapply(1:length(obs$Survey), function(x) ifelse(obs$Survey[x]=="Spawners", 1, ifelse(obs$Survey[x]=="Juveniles", 2, NA))))

Data_count_spawn <- Data_count %>% filter(Category == "Spawners")
Data_count_juv <- Data_count %>% filter(Category == "Juveniles")

Data_dens <- data.frame( "Catch_KG" = obs$Density, 
              "Year" = as.numeric(obs$Year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1,
               "Lat" = obs$Lat, 
               "Lon" = obs$Lon, 
               "Pass" = 0,
               "Knot" = unlist(obs$child_i),
               "Category" = obs$Survey,
               "CategoryNum" = sapply(1:length(obs$Survey), function(x) ifelse(obs$Survey[x]=="Spawners", 1, ifelse(obs$Survey[x]=="Juveniles", 2, NA))))

plot_network(Network_sz_LL = Network_sz_LL, Data = Data_count, arrows=TRUE, FilePath=fig_dir, FileName = "Survey_Locs")
plot_network(Network_sz_LL = Network_sz_LL, Data = Data_count, arrows=TRUE, byYear=TRUE, FilePath=fig_dir, FileName = "Survey_Locs_byYear")
plot_network(Network_sz_LL = Network_sz_LL, Data = Data_dens, arrows=TRUE, byYear=TRUE, byValue=TRUE, FilePath=fig_dir, FileName = "Density_byYear")
plot_network(Network_sz_LL = Network_sz_LL, Data = Data_dens %>% filter(Category=="Spawners"), value_label = "Spawner density", obs_color=RColorBrewer::brewer.pal(3,"Set1")[2], arrows=TRUE, byYear=TRUE, byValue=TRUE, FilePath=fig_dir, FileName = "Spawner_Density_byYear")
plot_network(Network_sz_LL = Network_sz_LL, Data = Data_dens %>% filter(Category=="Juveniles"), value_label = "Juvenile density", obs_color=RColorBrewer::brewer.pal(3,"Set1")[1], arrows=TRUE, byYear=TRUE, byValue=TRUE, FilePath=fig_dir, FileName = "Juvenile_Density_byYear")

####################
## Habitat
####################
habvar <- unique(hab$variable)
n_p <- length(habvar)
var_names <- c("Land cover", "Coho distribution", "Gradient", "Secondary channel area", "Large wood volume", "Large wood volume\nper 100m", "Primary channel area", "Pools per 100m", "Percent slackwater pools", "Weighted gravel in riffles")

## interpolation
hablist <- lapply(1:n_p, function(p){
	sub <- unique(hab %>% filter(variable == habvar[p])) #%>% select(-Year)) #%>% filter(Population=="Siletz")

	## if habitat variable is not land cover or coho distribution
	if(habvar[p] %in% c("land_cover", "Coho_distr") == FALSE){
		sub <- unique(sub %>% filter(Year > 0)) %>% select(-Year)
		sub$value <- as.numeric(sub$value)

		if(nrow(sub)<nrow(Network_sz)){
			## smooth across space

			interp_lat <- sub$Lat
			interp_lon <- sub$Lon
			interp_z <- sub$value				

			find_lat <- Network_sz_LL$Lat[which(Network_sz_LL$child_s %in% sub$child_i == FALSE)]
			find_lon <- Network_sz_LL$Lon[which(Network_sz_LL$child_s %in% sub$child_i == FALSE)]
			find_child <- Network_sz_LL$child_s[which(Network_sz_LL$child_s %in% sub$child_i == FALSE)]			

			# compute <- akima::interp(x = interp_lon, y = interp_lat, z = interp_z, xo = find_lon, yo=find_lat, linear=FALSE, extrap=TRUE, duplicate = "mean")
			compute <- akima::interpp(x = interp_lat, y = interp_lon, z = interp_z, xo=find_lat, yo=find_lon, duplicate = "mean", extrap=TRUE)			

			interp_df <- data.frame('Population'=unique(sub$Population), 'child_s'=find_child, 'Lon'=find_lon, 'Lat'=find_lat, 'HabitatImpact'=unique(sub$HabitatImpact), 'variable'=unique(sub$variable), 'value'=as.numeric(compute$z))				

			# p2 <- ggplot(interp_df)+
			# geom_point(aes(x=Lon,y=Lat,color=value)) +
			# geom_point(data=sub, aes(x=Lon, y=Lat, fill=value), pch=22, cex=3) +
			# guides(fill = FALSE) +
			# ggtitle(paste0(var_names[p], " (", habvar[p], ")", " interpolated (v1)")) +
			# xlab("Longitude") + ylab("Longitude") +
			# scale_color_viridis_c() +
			# scale_fill_viridis_c() +
			# mytheme()	
			# ggsave(file.path(fig_dir, paste0(habvar[p], "_interpolated_v1.png")),p2)		

			hab_info <- rbind.data.frame(interp_df, sub)
			obs_new <- sapply(1:nrow(hab_info), function(x) ifelse(is.na(hab_info$value[x]),mean(hab_info$value,na.rm=TRUE),hab_info$value[x]))
			hab_info$value <- obs_new			

			p3 <- ggplot(hab_info)+
			geom_point(aes(x=Lon,y=Lat,color=value)) +
			geom_point(data=sub, aes(x=Lon, y=Lat, fill=value), pch=22, cex=3) +
			guides(fill = FALSE, color = guide_legend(title = var_names[p])) +
			ggtitle(paste0(var_names[p], " (", habvar[p], ")", " interpolated")) +
			xlab("Longitude") + ylab("Longitude") +
			scale_color_viridis_c() +
			scale_fill_viridis_c() +
			mytheme()	
			ggsave(file.path(fig_dir, paste0(habvar[p], "_interpolated.png")),p3)		

			hab_new <- lapply(1:nrow(Network_sz_LL), function(x){
				child <- Network_sz_LL$child_s[x]
				find_hab <- hab_info %>% filter(child_s==child)
				if(nrow(find_hab)==1) return(find_hab)
				if(nrow(find_hab)>1){
					find_hab$value <- mean(find_hab$value)
					return(find_hab[1,])
				}
			})
			hab_new <- do.call(rbind, hab_new)

			sub <- hab_new			
			sub$value <- as.numeric(sub$value)

		}
	## if habitat variable is land cover or coho distribution
	} else {
		sub <- unique(sub %>% select(-Year))

		missing_child <- net_nodes[which(net_nodes %in% sub$child_s == FALSE)]
		if(length(missing_child)>0){
				fill <- lapply(1:length(missing_child), function(x){
					find_next <- Network_sz_LL %>% filter(child_s == missing_child[x])
					find_prev <- Network_sz_LL %>% filter(parent_s == missing_child[x])		

					if(nrow(find_next)>0){
						hab_next <- sub %>% filter(child_s == find_next$parent_s)
					} else { hab_next <- NULL }
					if(nrow(find_prev)>0){
						hab_prev <- sub %>% filter(child_s == find_prev$child_s)
					} else{ hab_prev <- NULL}		

					hab_fill <- NULL
					if(all(is.null(hab_next))==FALSE){
						if(nrow(hab_next)>0){
								hab_fill <- hab_next
								hab_fill$child_s = missing_child[x]
								hab_fill$Lat = Network_sz_LL[which(Network_sz_LL$child_s == missing_child[x]),"Lat"]
								hab_fill$Lon = Network_sz_LL[which(Network_sz_LL$child_s == missing_child[x]),"Lon"]
						}
					}
					if(all(is.null(hab_prev))==FALSE){
							if(nrow(hab_prev)>0){
								hab_fill <- hab_prev
								hab_fill$child_s = missing_child[x]
								hab_fill$Lat = Network_sz_LL[which(Network_sz_LL$child_s == missing_child[x]),"Lat"]
								hab_fill$Lon = Network_sz_LL[which(Network_sz_LL$child_s == missing_child[x]),"Lon"]
							}
					} 
					return(hab_fill)
				})	
				fill <- do.call(rbind, fill)

			p3 <- ggplot(fill)+
			geom_point(aes(x=Lon,y=Lat,color=value)) +
			geom_point(data=sub, aes(x=Lon, y=Lat, fill=value), pch=22, cex=3) +
			guides(color = FALSE, fill = guide_legend(title = var_names[p])) +
			ggtitle(paste0(var_names[p], " (", habvar[p], ")", " interpolated")) +
			xlab("Longitude") + ylab("Longitude") +
			scale_color_viridis_d() +
			scale_fill_viridis_d() +
			mytheme()	
			ggsave(file.path(fig_dir, paste0(habvar[p], "_interpolated.png")),p3)		
				
				sub <- rbind(sub, fill)
		}
	}
	return(sub)
})
hab_df <- do.call(rbind, hablist)

hab_sub <- hab_df %>% filter(variable == "land_cover")
land_cover_names <- unique(hab_sub$value)
land_cover_num <- seq_along(land_cover_names)
hab_sub$value <- sapply(1:nrow(hab_sub), function(x) land_cover_num[which(land_cover_names==hab_sub$value[x])])
hab_df[which(hab_df$variable == "land_cover"),"value"] <- hab_sub$value
hab_df <- unique(hab_df)

hab_spread <- hab_df %>% select(-c(Lon, Lat, HabitatImpact)) %>% tidyr::spread(key = variable, value = value)
hab_spread$GRADIENT <- as.numeric(hab_spread$GRADIENT)
hab_spread$LWDVOL1 <- as.numeric(hab_spread$LWDVOL1)
hab_spread$PCTSWPOOL <- as.numeric(hab_spread$PCTSWPOOL)
hab_spread$POOLS100 <- as.numeric(hab_spread$POOLS100)
hab_spread$PRICHNAREA <- as.numeric(hab_spread$PRICHNAREA)
hab_spread$SECCHNAREA <- as.numeric(hab_spread$SECCHNAREA)
hab_spread$VOLUMELWD <- as.numeric(hab_spread$VOLUMELWD)
hab_spread$WGTED_SLOPE_GRAVEL <- as.numeric(hab_spread$WGTED_SLOPE_GRAVEL)

require(GGally)
ppairs <- ggpairs(data=hab_spread, columns=3:ncol(hab_spread)) +
		mytheme()

hab_spread2 <- hab_spread %>% select(-c(PCTSWPOOL, LWDVOL1))
ppairs <- ggpairs(data=hab_spread2, columns=3:ncol(hab_spread2)) +
		mytheme()

habvar_df <- data.frame("variable"=habvar, 
						"HabitatImpact"=sapply(1:length(habvar), function(x) unique(unlist(hab_df %>% filter(variable==habvar[x]) %>% select(HabitatImpact)))))
habvar_df$toUse <- 1
habvar_df[which(habvar_df$variable %in% c("PCTSWPOOL", "LWDVOL1")),"toUse"] <- 0
habvar_df$HabitatImpact <- as.character(habvar_df$HabitatImpact)
habvar_df$HabitatImpact[which(habvar_df$HabitatImpact=="Coho")] <- "SpawnersJuveniles"

## Setup habitat covariates
## number of network nodes
n_x <- nrow(Network_sz_LL)

## number of years
n_t <- length(min(Data_count$Year):max(Data_count$Year))
n_t_spawn <- length(min(Data_count_spawn$Year):max(Data_count_spawn$Year))
n_t_juv <- length(min(Data_count_juv$Year):max(Data_count_juv$Year))

## number of habitat covariates
n_p <- length(habvar)

## number of observations
n_i <- nrow(Data_count)
n_i_spawn <- nrow(Data_count_spawn)
n_i_juv <- nrow(Data_count_juv)

### all data
X_gtp_all <- array(0, dim=c(n_x,n_t,n_p))
for(p in 1:n_p){
	psub <- hab_df %>% filter(variable == habvar[p])
	mat <- matrix(0, nrow=n_x, ncol=1)
	mat[psub$child_s,1] <- as.numeric(psub$value)
	mat_sd <- (mat - mean(mat))/sd(mat)
	X_gtp_all[,,p] <- mat_sd
}

X_itp_all <- array(0, dim=c(n_i,n_t,n_p))
for(i in 1:n_i){
	for(p in 1:n_p){
		knot <- Data_count$Knot[i]
		index <- which(net_nodes == knot)
		X_itp_all[i,,p] <- X_gtp_all[index,,p]
	}
}

Xconfig_all <- array(1, dim=c(2,2,n_p))
## remove related
Xconfig_all[,,which(habvar_df$toUse==0)] <- 0
## remove things that don't apply to spawners
Xconfig_all[,1,which(grepl("Spawners",habvar_df$HabitatImpact)==FALSE)] <- 0
## remove things that don't apply to juveniles
Xconfig_all[,2,which(grepl("Juveniles",habvar_df$HabitatImpact)==FALSE)] <- 0

## more basic habitat variables
Xconfig_all2 <- Xconfig_all
Xconfig_all2[,,which(habvar %in% c("land_cover", "Coho_distr")==FALSE)] <- 0


## spawners
X_gtp_spawn <- array(0, dim=c(n_x,n_t_spawn,n_p))
for(p in 1:n_p){
	psub <- hab_df %>% filter(variable == habvar[p])
	mat <- matrix(0, nrow=n_x, ncol=1)
	mat[psub$child_s,1] <- as.numeric(psub$value)
	mat_sd <- (mat - mean(mat))/sd(mat)
	X_gtp_spawn[,,p] <- mat_sd
}

X_itp_spawn <- array(0, dim=c(n_i_spawn,n_t_spawn,n_p))
for(i in 1:n_i_spawn){
	for(p in 1:n_p){
		knot <- Data_count_spawn$Knot[i]
		index <- which(net_nodes == knot)
		X_itp_spawn[i,,p] <- X_gtp_spawn[index,,p]
	}
}

Xconfig_spawn <- array(Xconfig_all[,1,], dim=c(dim(Xconfig_all)[1], 1, dim(Xconfig_all)[3]))
Xconfig_spawn2 <- Xconfig_spawn
Xconfig_spawn2[,,which(habvar %in% c("land_cover", "Coho_distr")==FALSE)] <- 0

## juveniles
X_gtp_juv <- array(0, dim=c(n_x,n_t_juv,n_p))
for(p in 1:n_p){
	psub <- hab_df %>% filter(variable == habvar[p])
	mat <- matrix(0, nrow=n_x, ncol=1)
	mat[psub$child_s,1] <- as.numeric(psub$value)
	mat_sd <- (mat - mean(mat))/sd(mat)
	X_gtp_juv[,,p] <- mat_sd
}

X_itp_juv <- array(0, dim=c(n_i_juv,n_t_juv,n_p))
for(i in 1:n_i_juv){
	for(p in 1:n_p){
		knot <- Data_count_juv$Knot[i]
		index <- which(net_nodes == knot)
		X_itp_juv[i,,p] <- X_gtp_juv[index,,p]
	}
}

Xconfig_juv <- array(Xconfig_all[,2,], dim=c(dim(Xconfig_all)[1], 1, dim(Xconfig_all)[3]))
Xconfig_juv2 <- Xconfig_juv
Xconfig_juv2[,,which(habvar %in% c("land_cover", "Coho_distr")==FALSE)] <- 0


########################
## spawners_habitat_discrete
#########################
path <- file.path(sil_dir, 'spawners_habitat_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count_spawn
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE,
                  X_gtp = X_gtp_spawn, X_itp = X_itp_spawn,
                  Xconfig_zcp = Xconfig_spawn)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))
Map[["gamma1_ctp"]] <- factor(rep(NA, length(Par[["gamma1_ctp"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_spawn, X_itp = X_itp_spawn, 
                  Xconfig_zcp = Xconfig_spawn,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_spawn, X_itp = X_itp_spawn,
                  Xconfig_zcp = Xconfig_spawn)

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names="Spawners", cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names="Spawners")
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names="Spawners")

########################
## spawners_habitat_discrete
#########################
path <- file.path(sil_dir, 'spawners_basichab_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count_spawn
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE,
                  X_gtp = X_gtp_spawn, X_itp = X_itp_spawn,
                  Xconfig_zcp = Xconfig_spawn2)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))
Map[["gamma1_ctp"]] <- factor(rep(NA, length(Par[["gamma1_ctp"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_spawn, X_itp = X_itp_spawn, 
                  Xconfig_zcp = Xconfig_spawn2,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_spawn, X_itp = X_itp_spawn,
                  Xconfig_zcp = Xconfig_spawn2)

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names="Spawners", cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names="Spawners")
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names="Spawners")


########################
## spawners_discrete
#########################
path <- file.path(sil_dir, 'spawners_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count_spawn
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map))

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names="Spawners", cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names="Spawners")
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names="Spawners")



########################
## juveniles_habitat_discrete
#########################
path <- file.path(sil_dir, 'juveniles_habitat_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count_juv
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=0)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE,
                  X_gtp = X_gtp_juv, X_itp = X_itp_juv,
                  Xconfig_zcp = Xconfig_juv)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))
Map[["gamma1_ctp"]] <- factor(rep(NA, length(Par[["gamma1_ctp"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_juv, X_itp = X_itp_juv,
                  Xconfig_zcp = Xconfig_juv,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_juv, X_itp = X_itp_juv,
                  Xconfig_zcp = Xconfig_juv)

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names="Spawners", cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names="Spawners")
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names="Spawners")

########################
## juveniles_basichab_discrete
#########################
path <- file.path(sil_dir, 'juveniles_basichab_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count_juv
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=0)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE,
                  X_gtp = X_gtp_juv, X_itp = X_itp_juv,
                  Xconfig_zcp = Xconfig_juv2)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))
Map[["gamma1_ctp"]] <- factor(rep(NA, length(Par[["gamma1_ctp"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_juv, X_itp = X_itp_juv,
                  Xconfig_zcp = Xconfig_juv2,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_juv, X_itp = X_itp_juv,
                  Xconfig_zcp = Xconfig_juv2)

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names="Spawners", cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names="Spawners")
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names="Spawners")

########################
## juveniles_discrete
#########################
path <- file.path(sil_dir, 'juveniles_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count_juv
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=0)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map))

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names="Spawners", cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names="Spawners")
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names="Spawners")

########################
## multivar_habitat_discrete
#########################
path <- file.path(sil_dir, 'multivar_habitat_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=2, "Epsilon2"=0)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE,
                  X_gtp = X_gtp_all, X_itp = X_itp_all,
                  Xconfig_zcp = Xconfig_all)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))
Map[["gamma1_ctp"]] <- factor(rep(NA, length(Par[["gamma1_ctp"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_all, X_itp = X_itp_all,
                  Xconfig_zcp = Xconfig_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_all, X_itp = X_itp_all,
                  Xconfig_zcp = Xconfig_all)

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names=c("Spawners","Juveniles"), cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names=c("Spawners","Juveniles"))
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=c("Spawners","Juveniles"))

########################
## multivar_discrete
#########################
path <- file.path(sil_dir, 'multivar_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=2, "Epsilon2"=0)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))
# Map[["gamma1_ctp"]] <- factor(rep(NA, length(Par[["gamma1_ctp"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map))

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names=c("Spawners","Juveniles"), cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names=c("Spawners","Juveniles"))
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=c("Spawners","Juveniles"))

########################
## multivar_basichab_discrete
#########################
path <- file.path(sil_dir, 'multivar_basichab_discrete')
dir.create(path, showWarnings=FALSE)
setwd(path)

fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.dll"), to = path)
ignore <- file.copy(from = file.path(sil_dir, "VAST_v8_2_0.o"), to = path)


## spawners only
Data <- Data_count
# Data$Catch_KG[which(Data$Catch_KG > 0)] <- log(Data$Catch_KG[which(Data$Catch_KG > 0)])
Xconfig_inp <- array(0, dim=c(2,2,dim(X_gtp_all)[3]))
Xconfig_inp[,,1:2] <- 1

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=2, "Epsilon2"=0)

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=3, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=11,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, bias.correct=FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE,
                  X_gtp = X_gtp_all, X_itp = X_itp_all,
                  Xconfig_zcp = Xconfig_all2)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
# Map[["beta1_ft"]] <- factor(rep(NA, length(Par[["beta1_ft"]])))
Map[["gamma1_ctp"]] <- factor(rep(NA, length(Par[["gamma1_ctp"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0),
                  X_gtp = X_gtp_all, X_itp = X_itp_all,
                  Xconfig_zcp = Xconfig_all2)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=as.numeric(Data[,"CategoryNum"])-1, 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])), 
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  X_gtp = X_gtp_all, X_itp = X_itp_all,
                  Xconfig_zcp = Xconfig_all2)

fit$parameter_estimates$diagnostics

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(plot_set=c(1:7,11,13:14), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, TmbData=fit$data_list, spatial_list=fit$spatial_list, DirName=fig, category_names=c("Spawners","Juveniles"), cex=0.5)
plot_biomass_index(TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=fig, category_names=c("Spawners","Juveniles"))
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)
plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=c("Spawners","Juveniles"))

 df <- data.frame("Model"=c("multivar_discrete",
 							"multivar_basichab_discrete",
 							"multivar_habitat_discrete"))
 df$AIC <- NA
 for(i in 1:nrow(df)){
 	res <- readRDS(file.path(sil_dir, df[i,"Model"], "Fit.rds"))
 	df$AIC[i] <- res$parameter_estimates$AIC
 }
 df$deltaAIC <- df$AIC - min(df$AIC)
 df[order(df$deltaAIC),]


