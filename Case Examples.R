#BST 245 Project Script

#Load in states and their corresponding centers of populations (latitude/longitude)
states <- read.csv("state_pop_centers_coordinates.txt")
cancer <- read.csv("state_cancer_incidence_2015.csv") #lung cancer - clusters due to smoking
#cancer <- read.csv("state_skin_cancer_incidence_2015.csv") #skin cancer
    #https://www.statecancerprofiles.cancer.gov/map/map.withimage.php?00&006&053&00&2&01&1&1&10&0#results
#cancer <- read.csv("CurrentAdultSmoking.csv") #smoking
      #https://www.cdc.gov/statesystem/cigaretteuseadult.html



states$STNAME #Includes Puerto Rico and District of Columbia
cancer$State #includes Puerto Rico and District of Columbia

#Sort by state/territory name
states <- states[order(states$STNAME),] #US states and their centers of population according to 2010 census
  #https://www.census.gov/geographies/reference-files/2010/geo/2010-centers-population.html
cancer<- cancer[order(cancer$State),] #Lung cancer incidence per 100k among men 65+, all races, in 2015
  #https://www.statecancerprofiles.cancer.gov/map/map.withimage.php?00&157&047&00&1&01&1&1&10&0#results
sum(states$STNAME == cancer$State) #52 indicates alignement of state by row in each data set

#Make a new column of lower case names for mapping later on. Remove Alaska, Hawaii, and Puerto Rico
cancer$region <- tolower(cancer$State)
states$region <- tolower(states$STNAME)
cancer <- cancer[which(!(cancer$region %in% c("alaska", "hawaii", "puerto rico", "guam"))),]
states <- states[which(!(states$region %in% c("alaska", "hawaii", "puerto rico", "guam"))),]

#Remove DC as well
#cancer <- cancer[which(!(cancer$region %in% c("district of columbia"))),]
#states <- states[which(!(states$region %in% c("district of columbia"))),]

#Function to calculate distance between points given latitude/longitude
#https://blog.exploratory.io/calculating-distances-between-two-geo-coded-locations-358e65fcafae
get_geo_distance = function(long1, lat1, long2, lat2, units = "miles") {
  loadNamespace("purrr")
  loadNamespace("geosphere")
  longlat1 = purrr::map2(long1, lat1, function(x,y) c(x,y))
  longlat2 = purrr::map2(long2, lat2, function(x,y) c(x,y))
  distance_list = purrr::map2(longlat1, longlat2, function(x,y) geosphere::distHaversine(x, y))
  distance_m = sapply(distance_list, function(col) { col[1] }) #Line adjusted from original code, see resposnes at url
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

#geo distance between alabama and california
get_geo_distance(long1=states[1,5], lat1=states[1,4], long2=states[4,5], lat2=states[4,4])

#Get distance matrix between states/territories
dist <- matrix(NA, nrow=nrow(states), ncol=nrow(states))
for(i in 1:nrow(states)){
  for(j in 1:nrow(states)){
    dist[i, j] <- get_geo_distance(long1=states[i,5], lat1=states[i,4], long2=states[j,5], lat2=states[j,4])
  }
}


#Transform the distance matrix into binary weights where w_(ij) = 1 if the distance < d and w_(ij) = 0 otherwise
getBinaryWeights <- function(dist, d){
  filteredDist <- matrix(0, nrow=nrow(dist), ncol=ncol(dist))
  filteredDist[which(dist < d, arr.ind=T)] <- 1
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
  #Total number of observations
  n <- length(x)
  
  #Weights for Getis-Ord Gi*
  filteredDist <- getBinaryWeights(dist, d)
  
  #Intermediate values for getting results
  yi1.star <- mean(x)
  yi2.star <- (1/n)*sum(x^2) - yi1.star^2
  
  #Matrix to store results
  results <- matrix(NA, nrow=length(x), ncol=5)
  colnames(results) <- c("Gi*", "Expectation", "Variance", "z-stat", "p-value")
    
  #Get the relevant statistic and values for each observed value
  for(i in 1:n){
    weights <- filteredDist[i,]
    gi.star <- sum(weights*x)/sum(x)
    wi.star <- sum(weights)
    expectation <- mean(weights)
    variance <- (wi.star*(n-wi.star)*yi2.star)/((n^2)*(n-1)*(yi1.star^2))
    z.star <- (gi.star - expectation)/sqrt(variance)
    
    #Get two-tailed p-value
    p.value <- 2*pnorm(-abs(z.star), mean=0, sd=1)
    
    #Store results for i
    results[i,] <- c(gi.star, expectation, variance, z.star, p.value)
  }
  
  return(as.data.frame(results))
}

#Should yield uniformly distributed p-values
#results <- getGiStar(x=rnorm(49), dist=dist, d=1000)

#Lung cancer data
#getisOrdResults <- getGiStar(x=cancer$Age.Adjusted.Incidence.Rate, dist=dist, d=500)
#cancer$State[which(getisOrdResults$`p-value` < 0.05)] #States at 5% significance level
#cancer$State[which(getisOrdResults$`p-value` < 0.001)] #States at 5% significance level with correction
#numInRange(dist, d=500)
#states$STNAME[which(numInRange(dist, 500) < 5)] #Those states with few neighbors
#states$STNAME[which(numInRange(dist, 500) > 15)] #Those with many neighbors


#Returns array of Local Moran's I_i statistics, expected values, variances, z-stat, and p-values (assuming approximate normality)
#for every point supplied. If x is an attribute with n observations, dist should be an n by n matrix of 
#distances between observations with the same ordering as x, d should be the filtering distance
#IdentityWeightNull = T sets w_(ii) = 0 for all i. If false, w_(ii) takes the usual weight.
getLocalMoranIi <- function(x, dist, d, identityWeightNull = T){
  #Total number of observations
  n <- length(x)
  
  #Weights for Local Moran's Ii
  filteredDist <- getBinaryWeights(dist, d)
  
  #If identityWeightNull=T, set w_(ii) = 0 for all i
  if(identityWeightNull){
    diag(filteredDist) <- rep(0, n)
  }
  
  #Standardize the observed values
  z <- (x - mean(x))/sd(x)
  
  #Intermediate values for getting results
  b <- ((1/n)*sum(z^4))/(((1/n)*sum(z^2))^2)
  
  #Matrix to store results
  results <- matrix(NA, nrow=length(x), ncol=5)
  colnames(results) <- c("Ii", "Expectation", "Variance", "z-stat", "p-value")
  
  #return(z)
  
  #Get the relevant statistic and values for each observed value
  for(i in 1:n){
    weights <- filteredDist[i,]
    Istat <- z[i]*sum(weights*z)
    #print(paste("i=", i, "z=", z, "  Istat=", Istat))
    #print(paste("i=",i, "z=",z))
    expectation <- (-1/(n-1))*sum(weights)
    sumForVar <- 0
    for(j in 1:n){
      for(k in 1:n){
        sumForVar <- sumForVar + weights[j]*weights[k]
      }
    }
    variance <- (1/(n-1))*((n-b)*sum(weights^2)) + (1/((n-1)*(n-2)))*((2*b-n)*sumForVar) - (1/((n-1)^2))*sum(weights^2)
    z.stat <- (Istat - expectation)/sqrt(variance)
    
    #Get two-tailed p-value
    p.value <- 2*pnorm(-abs(z.stat), mean=0, sd=1)
    
    #Store results for i
    results[i,] <- c(Istat, expectation, variance, z.stat, p.value)
  }
  
  return(as.data.frame(results))
}

#Should yield uniformly distributed p-values
#results <- getLocalMoranIi(x=rnorm(49), dist=dist, d=500)
#summary(results$`p-value`)

#Lung cancer data
#localMoranResults <- getLocalMoranIi(x=cancer$Age.Adjusted.Incidence.Rate, dist=dist, d=500)
#cancer$State[which(localMoranResults$`p-value` < 0.05)] #States at 5% significance level
#cancer$State[which(localMoranResults$`p-value` < 0.001)] #States at 5% significance level with correction
#numInRange(dist, d=500)
#states$STNAME[which(numInRange(dist, 500) < 5)] #Those states with few neighbors
#states$STNAME[which(numInRange(dist, 500) > 15)] #Those with many neighbors





#Getting table for presentation of significant state clusters
#alpha_new = alpha/n = 0.05/49 = 0.00102 (Bonferroni) or
#alpha_new = 1 - (1-\alpha)^(1/n) = 0.00105 (Sidak)
dist <- matrix(NA, nrow=nrow(states), ncol=nrow(states))
for(i in 1:nrow(states)){
  for(j in 1:nrow(states)){
    dist[i, j] <- get_geo_distance(long1=states[i,5], lat1=states[i,4], long2=states[j,5], lat2=states[j,4])
  }
}

#Get global Moran I
#https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/
library(ape)
Moran.I(cancer$Age.Adjusted.Incidence.Rate, dist, na.rm=T)
#Moran.I(cancer$Data_Value, dist)

#d=500
getisOrdResults <- getGiStar(x=cancer$Age.Adjusted.Incidence.Rate, dist=dist, d=500)
localMoranResults <- getLocalMoranIi(x=cancer$Age.Adjusted.Incidence.Rate, dist=dist, d=500)
#getisOrdResults <- getGiStar(x=cancer$Data_Value, dist=dist, d=300)
#localMoranResults <- getLocalMoranIi(x=cancer$Data_Value, dist=dist, d=300)

getisSig <- which(getisOrdResults$`p-value` < 0.00102) #States at 5% significance level with correction
moranSig <- which(localMoranResults$`p-value` < 0.00102) #States at 5% significance level with correction
sig <- sort(unique(c(getisSig, moranSig)))

combinedResults <- cbind(getisOrdResults[,c("Gi*", "z-stat", "p-value")], localMoranResults[,c("Ii", "z-stat", "p-value")])
combinedResults$Region <- cancer$State
combinedResults <- combinedResults[,c(7, 1:6)]
combinedResults[sig,] 

library(xtable)
latex <- xtable(combinedResults[sig,], align=c("l", "l", rep("c", 6)), digits=c(0, 0, 2, 2, 3, 2, 2, 3))
print(latex, include.rownames=FALSE)




############################Generating a Map#########################################

#Separate Attempt
#https://gist.github.com/cdesante/4252133


#Third Attempt (works!!!)
#https://stackoverflow.com/questions/12341281/add-points-to-choropleth-map-in-ggplot2
all_states <- map_data("state")
all_states
head(all_states)
#cancer$region <- tolower(cancer$State)
#states$region <- tolower(states$STNAME)
Total <- merge(all_states, cancer, by="region")
head(Total)

#Remove alaska and hawaii and puerto Rico - nevermind, don't show up anyway
Total <- Total[which(!(Total$region %in% c("alaska", "hawaii", "puerto rico"))),]
#cancer <- cancer[which(!(cancer$region %in% c("alaska", "hawaii", "puerto rico"))),]
#states <- states[which(!(states$region %in% c("alaska", "hawaii", "puerto rico"))),]


#Warning message that appears is because Alaska, Hawaii, and Puerto Rico don't appear on the map
base_map<-qplot(long, lat, data=Total, group=group, fill = Age.Adjusted.Incidence.Rate, geom="polygon", xlim=c(-130, -65), ylim=c(25, 50)) + scale_fill_continuous(low = "darkblue", high = "red", guide="colorbar")
base_map
base_map + 
  geom_point(aes(LONGITUDE, LATITUDE,fill = NULL,group = NULL), size = 1,data=states, color="yellow") + 
  borders("state", size=.5) + 
  labs(fill = "Lung Cancer Incidence \nper 100,000" ,title = "Incidence Rates for Lung Cancer, 2015 (All Races, Male, Ages 65+) with 2010 Centers of Population", x="", y="") +
  theme_grey(base_size=16)  + xlab("Longitude") + ylab("Latitude")


#Same plot as third attempt but also with a circle for approximate 500 mile radius
#https://stackoverflow.com/questions/34183049/plot-circle-with-a-certain-radius-around-point-on-a-map-in-ggplot2
library(ggmap)
data = data.frame(
  ID = as.numeric(c(1:nrow(states))),
  longitude = states$LONGITUDE,
  latitude = states$LATITUDE)

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
myCircles <- make_circles(data, radius=500*1.60934, nPoints=200)
##################################################################################


#Generate plot as before but with circles
base_map<-qplot(long, lat, data=Total, group=group, fill = Age.Adjusted.Incidence.Rate, geom="polygon", xlim=c(-130, -65), ylim=c(25, 50), xlab="Longitude", ylab="Latitude") + scale_fill_continuous(low = "darkblue", high = "red", guide="colorbar")
new_map <- base_map + 
  geom_point(aes(LONGITUDE, LATITUDE,fill = NULL,group = NULL), size = 1,data=states, color="yellow") + 
  borders("state", size=.5) + 
  labs(fill = "Lung Cancer Incidence \nper 100,000" ,title = "Incidence Rates for Lung Cancer, 2015 (All Races, Male, Ages 65+) with 2010 Centers of Population", x="", y="") #+
  ########### add circles
new_map + geom_point(data = myCircles, aes(lon, lat, group = NULL, fill=NULL), color = "black", alpha = 0.1) +   theme_grey(base_size=16)

#Circle just around Kansas
which(states$region == "kansas")
new_map + geom_point(data = myCircles[myCircles$ID ==15,], aes(lon, lat, group = NULL, fill=NULL), color = "black", alpha = 0.2)+ 
  theme_grey(base_size=16) + xlab("Longitude") + ylab("Latitude")

#States within 500 miles of Kansas (check that distance calculation matches with map) - seems to be in close alignment
numInRange(dist, 500)[15]
states$region[which(dist[15,] < 50)]

