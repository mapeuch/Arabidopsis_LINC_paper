##################### Parameters ###########################

path<-"INSERT PATH/"
file<-"wt_results.csv" #choose either of the results file
output_name<-"wt"
bin_number<-15

##################### Libraries ############################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(stringr)
library(ggtext)
library(xfun)
library(ggthemes)
library(tidyverse)
library(svglite)
library(plyr)
library(grid)
library(gridExtra)

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

##################### Data prep ############################

df<-read.csv2(paste(path, file, sep=""))

tdf<-as.data.frame(t(df))
colnames(tdf)<-c(1:nrow(df))

tdf <- rownames_to_column(tdf, "marker")

tdf <- tdf %>%
  separate(marker, into = c("chromosome", "submarker"), sep = "_")
tdf$submarker<-as.numeric(tdf$submarker)

tdf <- tdf %>%
  mutate(chrom_id = case_when(
    chromosome == "K1" ~ 1,
    chromosome == "K2" ~ 2,
    chromosome == "K3" ~ 3,
    chromosome == "K4" ~ 4,
    chromosome == "K5" ~ 5
  ))
dim(tdf)

###################### exploring markers ########################

for (k in 1:5){
  K <- paste("K", k, sep="")
  K<-get_lines(tdf, c(K), tdf$chromosome )
  vector<-data.frame(value = K$submarker)
  plot<-ggplot(vector, aes(x = value, y = "value")) +
    geom_boxplot() +
    geom_point(color="blue", size=6)+
    ggtitle(paste("Marker distribution along chromosome", k, sep=" "))
  result <- plot
  assign(paste("plot_k", k, sep=""), result)
}

plot_k1
plot_k2
plot_k3
plot_k4
plot_k5

################## localizing COs ######################

# separate into 5 dfs

df_list <- split(tdf, tdf$chrom_id)
tdf1 <- df_list[[1]]  # dataframe where 'chromosome' is 1
tdf2 <- df_list[[2]]  # dataframe where 'chromosome' is 2
tdf3 <- df_list[[3]]  # dataframe where 'chromosome' is 3
tdf4 <- df_list[[4]]  # dataframe where 'chromosome' is 4
tdf5 <- df_list[[5]]  # dataframe where 'chromosome' is 5

# store the 5 dfs as raw tdfs

tdf_r1 <- tdf1
tdf_r2 <- tdf2
tdf_r3 <- tdf3
tdf_r4 <- tdf4
tdf_r5 <- tdf5

# replace NAs by genotype of previous marker
# might've been smarter to use closests marker rather than previous but I think 
# the distribution is pretty even so should not make a difference

for (k in 1:5){
  tdfi <- get(paste("tdf", k, sep=""))
  
  for (p in 3:ncol(tdfi)){
    for (i in 1:nrow(tdfi)){
      if (is.na(tdfi[i, p])){
        tdfi[i, p]<-tdfi[i,p-1]
      }
    }
  }
  result <- tdfi
  assign(paste("tdf", k, sep=""), result)
}


# create and add intervals to CO_df for each chromosome

plants <- ncol(tdf)-3
plant_id <- paste("plant", 1:plants, sep = "")
CO_df_cols<-c("chromosome", "interval", plant_id)

for (k in 1:5){
  tdfi <- get(paste("tdf", k, sep=""))
  len<-nrow(tdfi)-1

  CO_df <- data.frame(chromosome = numeric(len), interval = numeric(len))
  CO_df <- cbind(CO_df, data.frame(replicate(plants, rep(0, len))))
  colnames(CO_df) <- CO_df_cols

  for (i in 1:len){
    if (tdfi$submarker[i]<tdfi$submarker[i+1]){
      CO_df$chromosome[i]=k
      CO_df$interval[i]=paste(tdfi$submarker[i], tdfi$submarker[i+1], sep="_")
    }
  }
  
  result <- CO_df
  assign(paste("CO_df", k, sep=""), result)
}

# add 1 in intervals where COs occur
for (k in 1:5){
  CO_dfi <- get(paste("CO_df", k, sep=""))
  tdfi <- get(paste("tdf", k, sep=""))
  intervals<-nrow(CO_dfi)

  for (p in 3:c(plants+2)){
    for (i in 1:intervals){
        if (tdfi[i,p] != tdfi[i+1,p]){
          CO_dfi[i,p] = 1
      }
    }
  }
  result <- CO_dfi
  assign(paste("CO_df", k, sep=""), result)
  assign(paste("CO_df_raw", k, sep=""), result)
}


################## CO frequencies ######################

# individual interval frequencies
for (k in 1:5){
  CO_dfi <- get(paste("CO_df", k, sep=""))

  CO_dfi$occurence <- rowSums(CO_dfi[, 3:ncol(CO_dfi)])
  CO_dfi$frequency <- CO_dfi$occurence/plants
  hist(CO_dfi$frequency, main = paste("Distribution of the observed frequency of COs
     \n in each interval for chromosome",k), xlab = "Frequencies")
  CO_dfi$mean_position <- sapply(strsplit(as.character(CO_dfi$interval), "_"), function(x) mean(as.numeric(x)))
  plot<-ggplot(CO_dfi, aes(x = mean_position, y = occurence)) +
    geom_bar(stat = "identity", fill="lightslateblue")+
    geom_line(stat = "identity", col="hotpink1", size=1.2) +                
    scale_x_continuous(breaks = seq(0,
                                    round_any(max(CO_dfi$mean_position), 2000000, f = ceiling),
                                    by =2000000))+
    theme(axis.text.x = element_text(angle = 30, hjust = 1))+
    ggtitle(paste("CO position along chromosome", k, sep=" "))


  result <- CO_dfi
  assign(paste("CO_df", k, sep=""), result)
  assign(paste("CO_df_for_nanc", k, sep=""), result)
  assign(paste("plot", k, sep=""), plot)

}

#write.csv(CO_df1, paste(path,"CO_df1_wt_for_matlab.csv"))
plot1
plot2
plot3
plot4
plot5

# expected double frequencies
for (k in 1:5){
  CO_dfi <- get(paste("CO_df", k, sep=""))

  CO_dfi<-data.frame( chromosome=CO_dfi$chromosome,
                     interval=CO_dfi$interval,
                     position=CO_dfi$mean_position,
                     occurence=CO_dfi$occurence,
                     frequency=CO_dfi$frequency)
  intervals<-CO_dfi$interval
  mati <- matrix(0, nrow = length(intervals), ncol = length(intervals),
                dimnames = list(intervals, intervals))
  for(i in 1:nrow(CO_dfi)) {
    for(j in 1:nrow(CO_dfi)) {
      mati[CO_dfi$interval[i], CO_dfi$interval[j]] <- CO_dfi$frequency[i] * CO_dfi$frequency[j]
    }
  }

  assign(paste("CO_df", k, sep=""), CO_dfi)
  assign(paste("mat", k, sep=""), mati)
}

#write.csv(mat1, paste(path,"mat1_wt.csv"))

# observed double frequencies

for (k in 1:5){
  CO_dfi <- get(paste("CO_df_raw", k, sep=""))
  mati <- get(paste("mat", k, sep=""))

  sum_dfi <- data.frame(CO_numbers=(colSums(CO_dfi[,3:c(plants+2)])))
  intervals<-CO_dfi$interval
  mat_o <- matrix(0, nrow = length(intervals), ncol = length(intervals),
                    dimnames = list(intervals, intervals))
  for (p in 1:plants){
    COs<-sum_dfi[p,1]
    if (COs>1){
      position_df<-get_lines(CO_dfi[,c(2,p+2)], c(1), CO_dfi[,p+2])
      positions<-position_df$interval
      vec <- c(1:COs)
      combinations <- expand.grid(vec, vec)
      combinations <- combinations[combinations$Var1 != combinations$Var2, ]
      for (i in 1:nrow(combinations)){
        mat_o[positions[combinations[i,1]],positions[combinations[i,2]]]= mat_o[positions[combinations[i,1]],positions[combinations[i,2]]]+1
      }
    }
  }
  mat_o_raw<-mat_o
  mat_o<-mat_o/plants
  assign(paste("mat_o", k, sep=""), mat_o)
  assign(paste("mat_o_raw", k, sep=""), mat_o_raw)

}

#write.csv(mat_o_raw1, paste(path,"mat_o1_raw_wt.csv"))

################## calculate CoC ######################

# this is done separately for each chromosome
for (k in 1:5){
  mat_o <- get(paste("mat_o", k, sep=""))
  mat <- get(paste("mat", k, sep=""))
  result <- mat_o / mat
  assign(paste("CoC_matrix", k, sep="_"), result)
}

#write.csv(CoC_matrix_1, paste(path,"CoC_matrix_1_wt.csv"))

################## present data ######################

# make into long dfs
for (k in 1:5){
  mat_o <- get(paste("CoC_matrix", k, sep="_"))
  df_o <- as.data.frame(mat_o)
  df_o$position2 <- rownames(df_o)
  long_df <- melt(df_o, id.vars = 'position2')
  colnames(long_df)<-c("position1", "position2","CoC")
  result <- long_df
  assign(paste("Long_CoC", k, sep="_"), result)
}

# transform intervals into positions
for (k in 1:5){
  CoC <- get(paste("Long_CoC", k, sep="_"))
  CoC$mean_position1 <- sapply(strsplit(as.character(CoC$position1), "_"), function(x) mean(as.numeric(x)))
  CoC$mean_position2 <- sapply(strsplit(as.character(CoC$position2), "_"), function(x) mean(as.numeric(x)))
  result <- CoC
  assign(paste("Long_CoC", k, sep="_"), result)
}

# calculate distance between the two positions
for (k in 1:5){
  CoC <- get(paste("Long_CoC", k, sep="_"))
  CoC$distance <- abs(CoC$mean_position1 - CoC$mean_position2)
  result <- CoC
  assign(paste("Long_CoC", k, sep="_"), result)
}

#write.csv(Long_CoC_1, paste(path,"Long_CoC_1_wt.csv"))

# plot per chromosome
for (k in 1:5){
  CoC <- get(paste("Long_CoC", k, sep="_"))
  plot<-ggplot(CoC, aes(x = distance, y = CoC)) +
    coord_cartesian(ylim=c(0,3))+
    geom_point()+
    geom_smooth()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1))+
    ggtitle(paste("Chromosome", k, sep=" "))
  result <- plot
  assign(paste("CoC_plot", k, sep="_"), result)
}

# make it into a panel of plots for CoC plot per chromosome
grid.arrange(CoC_plot_1, CoC_plot_2, CoC_plot_3, CoC_plot_4, CoC_plot_5,
             top=textGrob("pss1 CoCs"))

# plot with all chromosomes together
long_full<-rbind(Long_CoC_1, Long_CoC_2, Long_CoC_3, Long_CoC_4, Long_CoC_5)

plot<-ggplot(long_full, aes(x = distance, y = CoC)) +
  geom_point()+
  coord_cartesian(ylim=c(0,3))+
  geom_smooth(se=FALSE)+
  ggtitle(file)
plot

################## averaging CoC ###########################

long_full<-rbind(Long_CoC_1, Long_CoC_2, Long_CoC_3, Long_CoC_4, Long_CoC_5)

# assign each CoC value to a bin based on the distance between the two positions
long_full<-long_full %>% mutate (bin = cut (distance, bin_number))

# get the mean position of these bins
mean<-tapply(long_full$distance, long_full$bin, mean)
mean_dis<-data.frame(bin=row.names(mean),
                    mean_dis=mean)
long_full <-merge(long_full, mean_dis, by="bin")

# get the mean CoC for each bin
mean<-tapply(long_full$CoC, long_full$bin, mean, na.rm=TRUE)
mean_CoC<-data.frame(bin=row.names(mean),
                     mean_CoC=mean)
long_full <-merge(long_full, mean_CoC, by="bin")

# get the CoC SD per bin
sd<-tapply(long_full$CoC, long_full$bin, sd, na.rm=TRUE)
sd_CoC<-data.frame(bin=row.names(sd),
                     sd_CoC=sd)
long_full <-merge(long_full, sd_CoC, by="bin")

# group all the averaged data in a simplified df
mean_df<-data.frame(mean_dis=long_full$mean_dis,
                    mean_CoC=long_full$mean_CoC,
                    sd_CoC=long_full$sd_CoC)
mean_df<-mean_df[!duplicated(mean_df),]

plot<-ggplot(mean_df, aes(x = mean_dis/10**6, y = mean_CoC)) +
  geom_point()+
  coord_cartesian(ylim=c(0,3))+
  geom_smooth()+
  ylab("Mean CoC")+
  xlab("Inter-interval distance (Mb)")+
  ggtitle(file)
plot

# save plot
svg(file=paste(path, paste(output_name,"CoC.svg", sep="_"), sep="/"))
plot
dev.off()
