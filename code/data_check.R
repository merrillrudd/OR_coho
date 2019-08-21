rm(list=ls())

####################
## load packages
####################

devtools::install_github("merrillrudd/RuddR")

library(tidyverse)
library(RuddR)

###################
## Directories
###################

main_dir <- "C:\\merrill\\OR_coho"
data_dir <- file.path(main_dir, "data")

fig_dir <- file.path(main_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

###########################
## Read in data
###########################

spawn_raw <- read.csv(file.path(data_dir, "spawners_wNetwork.csv"))
juv_raw <- read.csv(file.path(data_dir, "juveniles_wNetwork.csv"))
hab_data_raw <- read.csv(file.path(data_dir, "habitat_wNetwork.csv"))
hab_siletz_raw <- read.csv(file.path(data_dir, "habitat_Siletz.csv"))
## read network from habitat data
net_fromHab <- unique(hab_data_raw %>% select(Population, ParentNode, ChildNode, Shape_Leng, Node_WGSX, Node_WGSY, Coho_distr) %>% rename(parent_s=ParentNode, child_s=ChildNode, dist_s=Shape_Leng, Lon=Node_WGSX, Lat=Node_WGSY))

length(which(net_fromHab$Lat==0))
find_zero <- net_fromHab %>% filter(Lat==0)

## for now, plot region without zeros
net_sub <- net_fromHab %>% filter(Lat!=0)

pnet <- ggplot(net_sub) +
		geom_point(aes(x = Lon, y = Lat, color = Population)) +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
ggsave(file.path(fig_dir, "Oregon_coast.png"), pnet)

## read observations
obs_spawn <- spawn_raw %>% 
          filter(SpawningYear > 0) %>% 
          select(Population, SpawningYear, CohoAdultAUC, Survey_Length__km_, sampled_X, sampled_Y, ChildNode) %>% 
          rename(Year=SpawningYear, Count=CohoAdultAUC, dist_i=Survey_Length__km_, Lon=sampled_X, Lat=sampled_Y, child_i=ChildNode) %>% 
          mutate(Survey="Spawners") %>% 
          mutate(Density = Count/dist_i)# %>% select(-c(Count,dist_i))

obs_juv <- juv_raw %>% 
        filter(Year_ > 0) %>% 
        select(Population, Year_, CohoPerKil, sampled_X, sampled_Y, ChildNode, CohoParr, Survey_Len) %>% 
        rename(Year=Year_, Density=CohoPerKil, Count=CohoParr, dist_i=Survey_Len, Lat=sampled_Y, Lon=sampled_X, child_i=ChildNode) %>% 
        mutate(Survey="Juveniles")

obs_coho <- rbind.data.frame(obs_spawn, obs_juv)

pobs <- ggplot() +
	geom_point(data = net_sub, aes(x = Lon, y = Lat), color = "gray") +
	geom_point(data = obs_coho, aes(x = Lon, y = Lat, color = Survey)) +
	scale_color_brewer(palette = "Set1") +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Survey_locs.png"), pobs)

pobs2 <- ggplot() +
	geom_point(data = net_sub, aes(x = Lon, y = Lat), color = "gray") +
	geom_point(data = obs_coho, aes(x = Lon, y = Lat, color = Survey)) +
	scale_color_brewer(palette = "Set1") +
	facet_wrap(.~Population, scale = "free") +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Survey_locs_byRegion.png"), pobs2, width = 10, height=8)


### habitat info
hab_all <- hab_data_raw %>% select(Population, ChildNode, Node_WGSX, Node_WGSY, YEAR_, land_cover, Coho_distr, GRADIENT, SECCHNAREA, VOLUMELWD, LWDVOL1, PRICHNAREA) %>% mutate(HabitatImpact = "Coho")
hab_all_df <- hab_all %>% tidyr::gather(key = "variable", value = "value", land_cover:PRICHNAREA) %>%
			  rename(child_s = ChildNode, Lon = Node_WGSX, Lat = Node_WGSY, Year = YEAR_)

hab_juv <- hab_data_raw %>% select(Population, ChildNode, Node_WGSX, Node_WGSY, YEAR_, POOLS100, PCTSWPOOL) %>% mutate(HabitatImpact = "Juveniles")
hab_juv_df <- hab_juv %>% tidyr::gather(key = "variable", value = "value", POOLS100:PCTSWPOOL)%>%
			  rename(child_s = ChildNode, Lon = Node_WGSX, Lat = Node_WGSY, Year = YEAR_)

hab_spawn <- hab_siletz_raw %>% mutate(Population="Siletz") %>% select(Population, ChildNode, lat, long, YEAR_, WGTED_SLOPE_GRAVEL) %>% mutate(HabitatImpact = "Spawners")
hab_spawn_df <- hab_spawn %>% tidyr::gather(key = "variable", value = "value", WGTED_SLOPE_GRAVEL) %>%
			  rename(child_s = ChildNode, Lon = long, Lat = lat, Year = YEAR_)

hab <- unique(rbind.data.frame(hab_all_df, hab_juv_df , hab_spawn_df) %>% filter(Lat != 0)) #%>% filter(Year > 0) %>% select(-Year))

variables <- unique(hab$variable)
var_names <- c("Land cover", "Coho distribution", "Gradient", "Secondary channel area", "Large wood volume", "Large wood volume\nper 100m", "Primary channel area", "Pools per 100m", "Percent slackwater pools", "Weighted gravel in riffles")
for(i in 1:length(variables)){
	sub <- hab %>% filter(variable == variables[i]) #%>% filter(Population=="Siletz")
	if(variables[i] %in% c("land_cover", "Coho_distr") == FALSE){
		sub <- sub %>% filter(Year > 0)
		sub$value <- as.numeric(sub$value)
	}
	p <- ggplot(sub) +
		geom_point(aes(x = Lon, y = Lat, color = value), alpha=0.5) +
		guides(color = guide_legend(title = var_names[i])) +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()	
	if(length(unique(sub$value))<10) p <- p + scale_color_viridis_d()
	if(length(unique(sub$value))>10) p <- p + scale_color_viridis_c()
	ggsave(file.path(fig_dir, paste0("Hab_", variables[i], ".png")), p)
}



saveRDS(obs_coho, file.path(data_dir, "observations.rds"))
saveRDS(net_sub, file.path(data_dir, "network.rds"))
saveRDS(hab, file.path(data_dir, "habitat.rds"))
