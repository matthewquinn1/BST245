#BST 245 Project Script

#Working directory
directory <- "D:/Matt Work/Harvard/Multivariate and Longitudinal Data/Project/Multiple Cause of Death Data"
setwd(directory)

#Load in counties and their corresponding centers of populations (latitude/longitude)
counties <- read.csv("county_pop_centers_coordinates.txt")
#Load in data on drug-induced deaths per 100,000 individuals
deaths <- read.table("Multiple Cause of Death 2017.txt", sep="\t", header=TRUE) 

#Data cleaning
##Remove notes column and last row of the deaths data which has totals in it
deaths <- deaths[1:(nrow(deaths)-1), -1]

##Replace "Unreliable", "suppressed" and 0s with NA in deaths data
deaths[deaths == "Unreliable"] <- NA
deaths[deaths == "Suppressed"] <- NA
deaths[deaths == 0] <- NA


##Reformat the county names in the counties file to match the deaths file
##This is necessary for merging the two data frames and because there are counties in different
  ##states that share the same name
counties$County <- NA
for(i in 1:nrow(counties)){
  countyTemp <- counties$COUNAME[i]
  stateTemp <- counties$STNAME[i]
  counties$County[i] <- toString(paste0(countyTemp, " County, ", state.abb[which(state.name == stateTemp)]))
}

#Find duplicated county names
which(duplicated(counties$County))

##For counties data set counties that are actually cities or incorrect, adjust the entry
  ##1217 is Baltimore city
  ##1598 is St. Louis city
  ##2917 is Bedford city, VA
  ##2926 is Fairfax city, VA
  ##2928 is Franklin city, VA
  ##2946 is Richmond city, VA
  ##2947 is Roanoke city, VA
indices <- c(1217, 1598, 2917, 2926, 2928, 2946, 2947)
for(i in indices){
  countyTemp <- counties$COUNAME[i]
  stateTemp <- counties$STNAME[i]
  counties$County[i] <- toString(paste0(countyTemp, " city, ", state.abb[which(state.name == stateTemp)]))
}

#Convert Louisiana Counties to Parishes
louisiana <- which(counties$STNAME == "Louisiana")
for(i in louisiana){
  countyTemp <- counties$COUNAME[i]
  stateTemp <- counties$STNAME[i]
  counties$County[i] <- toString(paste0(countyTemp, " Parish, ", state.abb[which(state.name == stateTemp)]))
}
rm(louisiana)

##Now reformat the deaths data frame to have separate state and county columns as well, which will
  ##be helpfulf for generating maps
deaths$region <- deaths$subregion <- NA
for(i in 1:nrow(deaths)){
  regionTemp <- strsplit(as.character(deaths$County[i]), ", ")[[1]][2]
  regionTemp <- state.name[which(state.abb == regionTemp)]
  if(length(regionTemp) == 0){
    deaths$region[i] <- NA
  }
  else{
    deaths$region[i] <- regionTemp
  }
  
  subregionTemp <- strsplit(as.character(deaths$County[i]), ",")[[1]][1]
  subregionTemp <- gsub(" County", "", subregionTemp)
  subregionTemp <- gsub(" Parish", "", subregionTemp)
  subregionTemp <- gsub("[.]", "", subregionTemp)
  if(length(subregionTemp) == 0){
    deaths$subregion[i] <- NA
  }
  else{
    deaths$subregion[i] <- subregionTemp
  }
}

## Remove Alaska, Hawaii, and Puerto Rico for purposes of mapping
deaths <- deaths[which(!(deaths$region %in% c("Alaska", "Hawaii", "Puerto Rico"))),]
counties <- counties[which(!(counties$STNAME %in% c("Alaska", "Hawaii", "Puerto Rico"))),]

#Merge deaths and counties data sets 
mergedDF <- merge(deaths, counties[, c("LATITUDE", "LONGITUDE", "County")], by=c("County"), all.x=T)

#################DEFINING FUNCTIONS AND FINDING DISTANCE MATRIX######################

#Function to calculate distance between points given latitude/longitude
#https://blog.exploratory.io/calculating-distances-between-two-geo-coded-locations-358e65fcafae
get_geo_distance = function(long1, lat1, long2, lat2, units = "miles") {
  loadNamespace("purrr")
  loadNamespace("geosphere")
  longlat1 = purrr::map2(long1, lat1, function(x,y) c(x,y))
  longlat2 = purrr::map2(long2, lat2, function(x,y) c(x,y))
  distance_list = purrr::map2(longlat1, longlat2, function(x,y) geosphere::distHaversine(x, y))
  distance_m = sapply(distance_list, function(col) { col[1] }) #Line adjusted from original code, see responses at url
  if (units == "km") {
    distance = distance_m / 1000.0;
  }
  else if (units == "miles") {
    distance = distance_m / 1609.344
  }
  else {
    distance = distance_m
    # This will return in meter as same way as distHaversine function. 
  }
  distance
}

#Parallelized distance finding, can take a while to run
library(doSNOW)
library(foreach)
cl<-makeCluster(6) #change to your number of CPU cores
registerDoSNOW(cl)

dist <- foreach(i=1:nrow(mergedDF), .combine=cbind) %dopar% {
    if(is.na(mergedDF$Crude.Rate[i])){
      return(rep(NA, nrow(mergedDF)))
    }
    else{
      sapply(1:nrow(mergedDF), FUN=function(index) get_geo_distance(long1=mergedDF[i,"LONGITUDE"], lat1=mergedDF[i,"LATITUDE"], long2=mergedDF[index,"LONGITUDE"], lat2=mergedDF[index,"LATITUDE"]))
    }
} 
stopCluster(cl)

#Write matrix to a file if desired
#library(MASS)
#write.matrix(dist, file = "distance_matrix", sep = " ")

#Transform the distance matrix into binary weights where w_(ij) = 1 if the distance < d and w_(ij) = 0 otherwise
getBinaryWeights <- function(dist, d){
  filteredDist <- matrix(NA, nrow=nrow(dist), ncol=ncol(dist))
  filteredDist[which(dist <= d, arr.ind=T)] <- 1
  filteredDist[which(dist > d, arr.ind=T)] <- 0
  return(filteredDist)
}

#Function that returns the number of observations within d of a given observation, including itself
numInRange <- function(dist, d){
  filteredWeights <- getBinaryWeights(dist, d)
  apply(filteredWeights, MARGIN=1, FUN=sum)
}

#Returns array of Gi* statistics, expected values, variances, z-stat, and p-values (assuming approximate normality)
#for every point supplied. If x is an attribute with n observations, dist should be an n by n matrix of 
#distances between observations with the same ordering as x, d should be the filtering distance
getGiStar <- function(x, dist, d){
  x <- as.numeric(x)
  
  #Total number of complete observations is all observations minus those with a missing crude mortality
  #rate or a missing location (e.g. the distance to itself can't be calculated)
  n <- length(x) - sum(is.na(x) | is.na(diag(dist)))
  
  #Weights for Getis-Ord Gi*
  filteredDist <- getBinaryWeights(dist, d)
  
  #Intermediate values for getting results
  yi1.star <- mean(x, na.rm=T)
  yi2.star <- (1/n)*sum(x^2, na.rm=T) - yi1.star^2
  
  #Matrix to store results
  results <- matrix(NA, nrow=length(x), ncol=5)
  colnames(results) <- c("Gi*", "Expectation", "Variance", "z-stat", "p-value")
  
  #Get the relevant statistic and values for each observed value
  for(i in 1:length(x)){
    if(is.na(x[i]) || is.na(dist[i,i])){
      results[i,] <- c(NA, NA, NA, NA, NA)
    }
    else{
      weights <- filteredDist[i,]
      gi.star <- sum(weights*x, na.rm=T)/sum(x, na.rm=T)
      wi.star <- sum(weights, na.rm=T)
      expectation <- mean(weights, na.rm=T)
      variance <- (wi.star*(n-wi.star)*yi2.star)/((n^2)*(n-1)*(yi1.star^2))
      z.star <- (gi.star - expectation)/sqrt(variance)
      
      #Get two-tailed p-value
      p.value <- 2*pnorm(-abs(z.star), mean=0, sd=1)
      
      #Store results for i
      results[i,] <- c(gi.star, expectation, variance, z.star, p.value)
    }
  }
  
  return(as.data.frame(results))
}

#Should yield uniformly distributed p-values
#results <- getGiStar(x=rnorm(49), dist=dist, d=1000)

#Returns array of Local Moran's I_i statistics, expected values, variances, z-stat, and p-values (assuming approximate normality)
#for every point supplied. If x is an attribute with n observations, dist should be an n by n matrix of 
#distances between observations with the same ordering as x, d should be the filtering distance
#IdentityWeightNull = T sets w_(ii) = 0 for all i. If false, w_(ii) takes the usual weight.
getLocalMoranIi <- function(x, dist, d, identityWeightNull = T){
  x <- as.numeric(x)
  
  #Total number of complete observations is all observations minus those with a missing crude mortality
  #rate or a missing location (e.g. the distance to itself can't be calculated)
  n <- length(x) - sum(is.na(x) | is.na(diag(dist)))
  
  #Weights for Local Moran's Ii
  filteredDist <- getBinaryWeights(dist, d)
  
  #If identityWeightNull=T, set w_(ii) = 0 for all i for which location isn't NA
  if(identityWeightNull){
    diag(filteredDist)[!is.na(diag(filteredDist))] <- 0
  }
  
  #Standardize the observed values
  z <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  
  #Intermediate values for getting results
  b <- ((1/n)*sum(z^4, na.rm=T))/(((1/n)*sum(z^2, na.rm=T))^2)
  
  #Matrix to store results
  results <- matrix(NA, nrow=length(x), ncol=5)
  colnames(results) <- c("Ii", "Expectation", "Variance", "z-stat", "p-value")
  
  #return(z)
  
  #Get the relevant statistic and values for each observed value
  for(i in 1:length(x)){
    if(is.na(x[i]) || is.na(dist[i,i])){
      results[i,] <- c(NA, NA, NA, NA, NA)
    }
    else{
      weights <- filteredDist[i,]
      Istat <- z[i]*sum(weights*z, na.rm=T)
      #print(paste("i=", i, "z=", z, "  Istat=", Istat))
      #print(paste("i=",i, "z=",z))
      expectation <- (-1/(n-1))*sum(weights, na.rm=T)
      sumForVar <- 0
      sumForVar <- sum(sapply(1:length(weights), FUN = function(x) sum(weights[x]*weights, na.rm=T)), na.rm=T)
      # for(j in 1:length(weights)){
      #   if(!is.na(weights[j])){
      #     for(k in 1:length(weights)){
      #       if(!is.na(weights[k])){
      #         sumForVar <- sumForVar + weights[j]*weights[k]
      #       }
      #     }
      #   }
      # }
      variance <- (1/(n-1))*((n-b)*sum(weights^2, na.rm=T)) + (1/((n-1)*(n-2)))*((2*b-n)*sumForVar) - (1/((n-1)^2))*sum(weights^2, na.rm=T)
      z.stat <- (Istat - expectation)/sqrt(variance)
      
      #Get two-tailed p-value
      p.value <- 2*pnorm(-abs(z.stat), mean=0, sd=1)
      
      #Store results for i
      results[i,] <- c(Istat, expectation, variance, z.stat, p.value)
    }
    #print(i)
  }
  
  return(as.data.frame(results))
}

#Should yield uniformly distributed p-values
#results <- getLocalMoranIi(x=rnorm(49), dist=dist, d=500)
#summary(results$`p-value`)


##############RUNNING THE TESTS####################################

#Choice of d is somewhat arbitrary without subject-matter knowledge
getisOrdResults <- getGiStar(x=as.numeric(mergedDF$Crude.Rate), dist=dist, d=75)
localMoranResults <- getLocalMoranIi(x=as.numeric(mergedDF$Crude.Rate), dist=dist, d=75)
#getisOrdResults <- getGiStar(x=deaths$Data_Value, dist=dist, d=300)
#localMoranResults <- getLocalMoranIi(x=deaths$Data_Value, dist=dist, d=300)

cor(getisOrdResults$`p-value`, localMoranResults$`p-value`, use="pairwise.complete.obs") #About 0.67 correlation

#True n = nrow(mergedDF) - sum(is.na(mergedDF$Crude.Rate) | is.na(mergedDF$LATITUDE)) = 2280
numComplete <-  nrow(mergedDF) - sum(is.na(mergedDF$Crude.Rate) | is.na(mergedDF$LATITUDE)) 
#New alpha is:
bon <- 0.05/numComplete #0.0000219
sidak <- 1 - (1-0.05)^(1/numComplete) #0.0000225

getisSig <- which(getisOrdResults$`p-value` < sidak) #counties at 5% significance level with correction
moranSig <- which(localMoranResults$`p-value` < sidak) #counties at 5% significance level with correction
sig <- sort(unique(c(getisSig, moranSig)))
mergedDF[sig, "County"]
table(mergedDF[sig, "region"])

#Median p-value of Gi^* results that were not significant but are significant under I_i
#notSig <- moranSig[-c( which(getisSig %in% moranSig))]
#median(getisOrdResults$`p-value`[notSig])

combinedResults <- cbind(getisOrdResults[,c("Gi*", "z-stat", "p-value")], localMoranResults[,c("Ii", "z-stat", "p-value")])
combinedResults$Region <- mergedDF$County
combinedResults <- combinedResults[,c(7, 1:6)]
combinedResults[sig,] 

library(xtable)
latex <- xtable(combinedResults[sig,], align=c("l", "l", rep("c", 6)), digits=c(0, 0, 2, 2, 5, 2, 2, 5))
print(latex, include.rownames=FALSE)




############################Generating Maps#########################################
library(ggplot2)
library(mapproj)

#https://stackoverflow.com/questions/12341281/add-points-to-choropleth-map-in-ggplot2
all_counties <- map_data("county")
head(all_counties)

#Adjust deaths data set 
deaths$region <- tolower(deaths$region)
deaths$subregion <- tolower(deaths$subregion)

######ATTEMPT FOR COUNTIES - SEEMS TO WORK##########
#Modified from comment from Paulo at:
#https://stackoverflow.com/questions/23714052/ggplot-mapping-us-counties-problems-with-visualization-shapes-in-r
#Solution in comment is actually incorrect so it was modified here
m.usa <- map_data("county")
for(i in 1:nrow(m.usa)){
  countyTemp <- m.usa$subregion[i]
  stateTemp <- m.usa$region[i]
  m.usa$subregion[i] <- toString(paste0(countyTemp, ", ", stateTemp))
}
m.usa$id <- m.usa$subregion

#Make new df just to hold map fill values
df <- data.frame(region = unique(m.usa$subregion),
                 mortalityRate = rnorm(length(unique(m.usa$subregion)), 50, 10),
                 stringsAsFactors = F)

#Adjust deaths data set 
deaths$region <- tolower(deaths$region)
deaths$subregion <- tolower(deaths$subregion)

for(i in 1:nrow(deaths)){
  countyTemp <- deaths$subregion[i]
  stateTemp <- deaths$region[i]
  deaths$subregion[i] <- toString(paste0(countyTemp, ", ", stateTemp))
}

#Match up mortality rates from death data set to the county map data
for(i in 1:nrow(df)){
  temp <- deaths$Crude.Rate[which(deaths$subregion == df$region[i])]
  if(length(temp) == 0){
    df$mortalityRate[i] <- NA
  }
  else{
    df$mortalityRate[i] <- temp
  }
}

m.usa <- m.usa[ ,-5]
names(m.usa)[5] <- 'region'

#ggplot, for whatever reason, breaks if counties data doesn't have a region column
counties$region <- NA

#head(df)
#Plot of US counties with mortality rate as fill
base_plot <- ggplot(df, aes(map_id = region)) +
  geom_map(aes(fill = mortalityRate), map = m.usa) + 
  expand_limits(x = m.usa$long, y = m.usa$lat) +
  coord_map() + borders("state", colour="black") +
  labs(fill = "Drug-related Deaths \nper 100,000" ,title = "Crude Drug-related Mortality Rates, 2017 (All Races/Sexes/Ages) by County", x="", y="") +
  theme_grey(base_size=24)  + xlab("Longitude") + ylab("Latitude") + 
  scale_fill_gradient(low = "blue3", high = "red1")

png(filename="base_plot.png", width=1600, height=800)
base_plot
dev.off()


#Map with yellow dots indicating centers of population for each county
base_plot <- ggplot(df, aes(map_id = region)) +
  geom_map(aes(fill = mortalityRate), map = m.usa) + 
  expand_limits(x = m.usa$long, y = m.usa$lat) +
  coord_map() + borders("state", colour="black") +
  labs(fill = "Drug-related Deaths \nper 100,000" ,title = "Crude Drug-related Mortality Rates, 2017 (All Races/Sexes/Ages) with Centers of Population", x="", y="") +
  theme_grey(base_size=16)  + xlab("Longitude") + ylab("Latitude") + 
  scale_fill_gradient(low = "blue3", high = "red1") 
base_plot + geom_point(aes(LONGITUDE, LATITUDE), size = 0.5,data=counties, color="yellow")
  
png(filename="base_plot_with_centers.png", width=1600, height=800)
base_plot + geom_point(aes(LONGITUDE, LATITUDE), size = 0.5,data=counties, color="yellow")
dev.off()

##############################

#Same plot as before but also with a circle for approximate r mile radius
#https://stackoverflow.com/questions/34183049/plot-circle-with-a-certain-radius-around-point-on-a-map-in-ggplot2
library(ggmap)
data = data.frame(
  ID = as.numeric(c(1:nrow(counties))),
  longitude = counties$LONGITUDE,
  latitude = counties$LATITUDE)

#################################################################################
# create circles data frame from the centers data frame
make_circles <- function(centers, radius, nPoints = 100){
  # centers: the data frame of centers with ID
  # radius: radius measured in kilometer
  #
  meanLat <- mean(centers$latitude)
  # length per longitude changes with lattitude, so need correction
  radiusLon <- radius /111 / cos(meanLat/57.3) 
  radiusLat <- radius / 111
  circleDF <- data.frame(ID = rep(centers$ID, each = nPoints))
  angle <- seq(0,2*pi,length.out = nPoints)
  
  circleDF$lon <- unlist(lapply(centers$longitude, function(x) x + radiusLon * cos(angle)))
  circleDF$lat <- unlist(lapply(centers$latitude, function(x) x + radiusLat * sin(angle)))
  return(circleDF)
}

#Get the data frame of circles
#radius is in km, so do miles*1.60934
myCircles <- make_circles(data, radius=75*1.60934, nPoints=200)
##################################################################################


#Generate plot as before but with circles
base_map <- ggplot(df, aes(map_id = region)) +
  geom_map(aes(fill = mortalityRate), map = m.usa) + 
  expand_limits(x = m.usa$long, y = m.usa$lat) +
  coord_map() + borders("state", colour="black") +
  labs(fill = "Drug-related Deaths \nper 100,000" ,title = "Crude Drug-related Mortality Rates, 2017 (All Races/Sexes/Ages) with Centers of Population", x="", y="") +
  theme_grey(base_size=16)  + xlab("Longitude") + ylab("Latitude") + 
  scale_fill_gradient(low = "blue3", high = "red1") +
  geom_point(aes(LONGITUDE, LATITUDE), size = 1, data=counties, color="yellow")

#Circle just around Bronx County, NY
which(counties$COUNAME == "Bronx")
#Again need to add a region column to myCircles to avoid ggplot breaking
myCircles$region <- NA
base_map + geom_point(data = myCircles[myCircles$ID ==1797,], aes(lon, lat, group = NULL, fill=NULL), color = "green", alpha = 0.4, size=0.5)+ 
  theme_grey(base_size=16) + xlab("Longitude") + ylab("Latitude")

png(filename="base_plot_with_centers.png", width=1700, height=800)
base_map + geom_point(data = myCircles[myCircles$ID ==1797,], aes(lon, lat, group = NULL, fill=NULL), color = "green", alpha = 0.4, size=1.5)+ 
  theme_grey(base_size=24) + xlab("Longitude") + ylab("Latitude")
dev.off()



#######PLOTS INDICATING THOSE COUNTIES FOUND TO HAVE SIGNIFICANT SPATIAL CLUSTERING#######
library(ggplot2)
#Need the mergedDF combining the deaths and census data
head(mergedDF)
mergedDFTemp <- mergedDF

#Need sig to be defined from running the tests - see around line 290
head(sig)

#Append to the mergedDFTemp data set an indicator of statistical significance under at least one test
mergedDFTemp$significant <- NA
mergedDFTemp$significant[!is.na(getisOrdResults$`p-value`)] <- "No"
mergedDFTemp$significant[sig] <- "Yes"

m.usa <- map_data("county")
for(i in 1:nrow(m.usa)){
  countyTemp <- m.usa$subregion[i]
  stateTemp <- m.usa$region[i]
  m.usa$subregion[i] <- toString(paste0(countyTemp, ", ", stateTemp))
}
m.usa$id <- m.usa$subregion

#Make new df just to hold map fill values
df <- data.frame(region = unique(m.usa$subregion),
                 significant = rnorm(length(unique(m.usa$subregion)), 50, 10),
                 stringsAsFactors = F)

#Adjust mergedDFTemp data set 
mergedDFTemp$region <- tolower(mergedDFTemp$region)
mergedDFTemp$subregion <- tolower(mergedDFTemp$subregion)

for(i in 1:nrow(mergedDFTemp)){
  countyTemp <- mergedDFTemp$subregion[i]
  stateTemp <- mergedDFTemp$region[i]
  mergedDFTemp$subregion[i] <- toString(paste0(countyTemp, ", ", stateTemp))
}

#Match up significance from death data set to the county map data
for(i in 1:nrow(df)){
  temp <- mergedDFTemp$significant[which(mergedDFTemp$subregion == df$region[i])]
  if(length(temp) == 0){
    df$significant[i] <- NA
  }
  else{
    df$significant[i] <- temp
  }
}

m.usa <- m.usa[ ,-5]
names(m.usa)[5] <- 'region'

#ggplot, for whatever reason, breaks if counties data doesn't have a region column
counties$region <- NA

#head(df)
#Plot of US counties with mortality rate as fill
sig_plot <- ggplot(df, aes(map_id = region)) +
  geom_map(aes(fill = as.factor(significant)), map = m.usa) + 
  expand_limits(x = m.usa$long, y = m.usa$lat) +
  coord_map() + borders("state", colour="black") +
  labs(fill = "Significant" ,title = "Counties with Significant Crude Drug-related Mortality Rate Spatial Clustering, 2017", x="", y="") +
  theme_grey(base_size=24)  + xlab("Longitude") + ylab("Latitude") + 
  scale_fill_manual(values = c("thistle2", "red"), na.value="gray")

png(filename="significant_2017.png", width=1600, height=800)
sig_plot
dev.off()

