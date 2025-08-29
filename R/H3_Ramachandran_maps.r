#27/05/2021   Ivan Martin Hernandez



#Libraries
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())
library(dplyr)
library(stringr)
library(reshape2)
library(RColorBrewer)





#Import the data and make de DF
aa92 <- read.table(file = "ALLN92bin72.txt")
aa93 <- read.table(file = "ALLN93bin72.txt")
aa94 <- read.table(file = "ALLN94bin72.txt")
aa95 <- read.table(file = "ALLN95bin72.txt")
aa100 <- read.table(file = "ALLC100bin72.txt")
aa101 <- read.table(file = "ALLC101bin72.txt")
aa102 <- read.table(file = "ALLC102bin72.txt")
aa103 <- read.table(file = "ALLC103bin72.txt")




# Heat mapS


mid <- mean(aa92$V6)

a<-ggplot(aa92, aes(x = V2, y = V4)) + 
  geom_tile(aes(fill = V5),show.legend = FALSE)+
  #scale_fill_gradient(low = "white", high = "blue")+
  scale_fill_gradientn(colours = terrain.colors(7))+
  theme_linedraw()+
  theme(axis.ticks = element_blank(), axis.title = element_blank(),axis.text=element_text(size=14,face="bold"),
        panel.grid = element_blank())
print(a)




mid <- mean(aa103$V6)
b<-ggplot(aa103, aes(x = V2, y = V4)) + 
  geom_tile(aes(fill = V5),show.legend = FALSE)+
  #scale_fill_gradient2(midpoint = mid,low = "blue", mid = "white",high = "red")+
  scale_fill_distiller(palette = "Greys", direction = 1)+
  #scale_fill_distiller(palette = "RdGy", direction = -1)+
  theme_linedraw()+
  theme(axis.ticks = element_blank(), axis.title = element_blank(),axis.text=element_text(size=14,face="bold"),
        panel.grid = element_blank())
print(b)

