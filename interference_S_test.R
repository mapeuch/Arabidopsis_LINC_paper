##################### Parameters ###########################

# R script written by M.Puch perforiming test described in Jahns 2014 
#"Crossover Localisation Is Regulated by the Neddylation Posttranslational Regulatory Pathway"
# 10.1371/journal.pbio.1001930

sims <- 10000 # number of runs
path<-"INSERT PATH/"
file<-"sine3_MLH1_foci_positions.csv"

##################### Libraries ############################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)
library(dunn.test)
library(ggtext)
library(xfun)
library(ggthemes)
library(tidyverse)
library(svglite)
library(multcompView)
library(plyr)
library(grid)
library(gridExtra)
library(zoo)


get_lines <-function(df, vector, df_column){
  len<-nrow(df)
  result_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  colnames(result_df) <- colnames(df)
  for (i in 1:len){
    if (any(vector  %in% df_column[i])==TRUE){
      result_df<-rbind(result_df,df[i,])
    }
  }
  print(result_df)
}

##################### calculate S_exp ############################

df<-read.csv2(paste(path, file, sep=""))

df$SC_ID_tmp<-as.factor(paste(df$Cell_nb, df$SC_fragment, sep="_"))
df$SC_ID<-as.factor(paste(df$Slide, df$SC_ID_tmp, sep="_"))
df<-data.frame(SC_ID=df$SC_ID, 
               SC_length=df$SC_length, 
               MLH1_foci_nb=df$MLH1_foci_nb, 
               MLH1_foci_position=df$MLH1_foci_position)
fragments<-levels(df$SC_ID)

# get inter-CO distances

distance_df<-data.frame(SC_ID=c(),
                        Point1=c(),
                        Point2=c(),
                        Distance=c(),
                        MLH1_foci=c())

for (fragment in 1:length(fragments)){
  temp_df<-get_lines(df, c(fragments[fragment]), df$SC_ID)
  if (nrow(temp_df)>1){
      
    distances<-dist(temp_df$MLH1_foci_position)
    distance_matrix<-as.matrix(distances)
    
    tmp_dist_df <- as.data.frame(as.table(distance_matrix))
    colnames(tmp_dist_df) <- c("Point1", "Point2", "Distance")
    tmp_dist_df$Point1<-as.numeric(tmp_dist_df$Point1)
    tmp_dist_df$Point2<-as.numeric(tmp_dist_df$Point2)
    
    # reomve distances counted twice
    tmp_dist_df <- tmp_dist_df[tmp_dist_df$Point1 < tmp_dist_df$Point2, ]
    
    # keep only adjacent distances
    tmp_dist_df <- tmp_dist_df[tmp_dist_df$Point2 - tmp_dist_df$Point1 == 1, ]
    
    # add to general df
    tmp_dist_df_2<-data.frame(SC_ID=rep(temp_df[1,1], nrow(tmp_dist_df)),
                              Point1=tmp_dist_df$Point1,
                              Point2=tmp_dist_df$Point2,
                              Distance=tmp_dist_df$Distance,
                              MLH1_foci=rep(nrow(temp_df), nrow(tmp_dist_df)))
    distance_df<-rbind(distance_df, tmp_dist_df_2)
  }
}

hist(distance_df$Distance, breaks=20, xlab="adjacent inter-CO distances")
plot(density(distance_df$Distance), xlab="adjacent inter-CO distances")

# normalize these distances
norm_df<-data.frame(SC_ID=c(),
                    Point1=c(),
                    Point2=c(),
                    Distance=c(),
                    MLH1_foci=c(),
                    Normalized_distance=c())

for (fragment in 1:length(fragments)){
  
  temp_df<-get_lines(distance_df, c(fragments[fragment]), distance_df$SC_ID)
  
  if (nrow(temp_df)>1){
    cum_dist<-sum(temp_df$Distance)
    temp_df$Normalized_distance<-temp_df$Distance/cum_dist
    
    if (sum(temp_df$Normaized_distance !=1)){
      cat("Error: Element", fragment, "does not have a sum of normalized values equal to one.\n")
    }
    norm_df<-rbind(norm_df, temp_df)
  }
}

# calculate S exp
distance_df<-data.frame(SC_ID=norm_df$SC_ID,
                    Point1=norm_df$Point1,
                    Point2=norm_df$Point2,
                    MLH1_foci=norm_df$MLH1_foci,
                    Distance=norm_df$Normalized_distance)

fragments<-unique(distance_df$SC_ID)
k_vector<-rep(0, length(fragments))
S_vector<-rep(0, length(fragments))

for (fragment in 1:length(fragments)){
  
  temp_df<-get_lines(distance_df, c(fragments[fragment]), distance_df$SC_ID)
  k_vector[fragment]<-mean(temp_df$MLH1_foci)
  
  if (nrow(temp_df)>1){
    
    d=0
    for (i in 1:nrow(temp_df)) {
      d=d+(temp_df$Distance[i]-(1/nrow(temp_df)))**2
    }
    S_vector[fragment]=d*((nrow(temp_df)**2)*(nrow(temp_df)+1))/(nrow(temp_df)-1)
    
  }
}

S_exp=sum(S_vector)

##################### Simulate S_sim ############################

# Actual simulation part
k_df<-as.data.frame(table(k_vector)) # k_vector is the number of COs+1 in patches with 2 or more COs

simulation_df<-data.frame()

for (sim in 1:sims){
  pb <- txtProgressBar(min = 0, max = sims, style = 2)
  
  id=0
  for (k in k_vector){
    id=id+1
    random_decimals <- runif(k, min = 0, max = 1)
    df<-data.frame(Simulation = rep(sim, length(random_decimals)), ID=paste(id,k, sep="_"), Foci=c(1:length(random_decimals)),  Position = random_decimals)
    simulation_df <- rbind(simulation_df, df)
  }
  setTxtProgressBar(pb, sim)
}

# calculate the distance between the simulated COs

distance_df<-data.frame(Simulation=c(),
                        krun=c(),
                        Point1=c(),
                        Point2=c(),
                        Distance=c(),
                        foci=c())

for (sim in 1:sims){
  pb <- txtProgressBar(min = 0, max = sims, style = 2)
  
  temp_df_top<-get_lines(simulation_df, sim, simulation_df$Simulation)
  fragments<-unique(temp_df_top$ID)
  
for (fragment in 1:length(fragments)){
  temp_df<-get_lines(temp_df_top, c(fragments[fragment]), temp_df_top$ID)

    distances<-dist(temp_df$Position)
    distance_matrix<-as.matrix(distances)
    
    tmp_dist_df <- as.data.frame(as.table(distance_matrix))
    colnames(tmp_dist_df) <- c("Point1", "Point2", "Distance")
    tmp_dist_df$Point1<-as.numeric(tmp_dist_df$Point1)
    tmp_dist_df$Point2<-as.numeric(tmp_dist_df$Point2)
    
    # reomve distances counted twice
    tmp_dist_df <- tmp_dist_df[tmp_dist_df$Point1 < tmp_dist_df$Point2, ]
    
    # keep only adjacent distances
    tmp_dist_df <- tmp_dist_df[tmp_dist_df$Point2 - tmp_dist_df$Point1 == 1, ]
    
    # add to general df
    tmp_dist_df_2<-data.frame(Simulation=rep(sim, nrow(tmp_dist_df)),
                              krun=rep(temp_df[1,2], nrow(tmp_dist_df)),
                              Point1=tmp_dist_df$Point1,
                              Point2=tmp_dist_df$Point2,
                              Distance=tmp_dist_df$Distance,
                              foci=rep(nrow(temp_df), nrow(tmp_dist_df)))
    distance_df<-rbind(distance_df, tmp_dist_df_2)
  }
  setTxtProgressBar(pb, sim)
  
}

hist(distance_df$Distance, breaks=20, probability=TRUE)
plot(density(distance_df$Distance))

# normalize these distances
norm_df<-data.frame(Simulation=c(),
                    krun=c(),
                    Point1=c(),
                    Point2=c(),
                    Distance=c(),
                    Normalized_distance=c())

for (sim in 1:sims){
  pb <- txtProgressBar(min = 0, max = sims, style = 2)
  
  temp_df_top<-get_lines(distance_df, sim, distance_df$Simulation)
  fragments<-unique(temp_df_top$krun)
  
for (fragment in 1:length(fragments)){
  
  temp_df<-get_lines(temp_df_top, c(fragments[fragment]), temp_df_top$krun)
  
    cum_dist<-sum(temp_df$Distance)
    temp_df$Normalized_distance<-temp_df$Distance/cum_dist
    
    if (sum(temp_df$Normaized_distance !=1)){
      cat("Error: Element", fragment, "does not have a sum of normalized values equal to one.\n")
    }
    norm_df<-rbind(norm_df, temp_df)
  }
  setTxtProgressBar(pb, sim)
}

# calculate S_sim

distance_df<-data.frame(Simulation=norm_df$Simulation,
                        SC_ID=norm_df$krun,
                        Point1=norm_df$Point1,
                        Point2=norm_df$Point2,
                        MLH1_foci=norm_df$foci,
                        Distance=norm_df$Normalized_distance)

S_sim_vector<-rep(0, sims)

for (sim in 1:sims){
  pb <- txtProgressBar(min = 0, max = sims, style = 2)
  
  temp_df_top<-get_lines(distance_df, sim, distance_df$Simulation)
  fragments<-unique(distance_df$SC_ID)
  S_vector<-rep(0, length(fragments))

for (fragment in 1:length(fragments)){
  
  temp_df<-get_lines(temp_df_top, c(fragments[fragment]), temp_df_top$SC_ID)

    d=0
    for (i in 1:nrow(temp_df)) {
      d=d+(temp_df$Distance[i]-(1/nrow(temp_df)))**2
    }
    S_vector[fragment]=d*((nrow(temp_df)**2)*(nrow(temp_df)+1))/(nrow(temp_df)-1)
    
}
  S_sim_vector[sim]=sum(S_vector)
  setTxtProgressBar(pb, sim)
}

S_sims_below_S_exp<-0

for (S in S_sim_vector){
  if (S<S_exp){
    S_sims_below_S_exp=S_sims_below_S_exp+1
  }
}

p_value=S_sims_below_S_exp/sims

# save simulation data

setwd(path)

write.csv(simulation_df, paste(sim,"Simulation_data.csv", sep="_"), row.names = FALSE)
write.csv(distance_df, paste(sim,"Simulation_data.csv", sep="_"), row.names = FALSE)
write.csv(S_sim_vector, paste(sim,"S_sim_values.csv", sep="_"), row.names = FALSE)


