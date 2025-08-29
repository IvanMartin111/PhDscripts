#20/11/2022   Ivan Martin Hernandez
#Script that paints Ramachandram maps from files of Phi/Phi anles


#Libraries
library(ggplot2)
library(ggrepel)
library(ggpubr)
theme_set(theme_pubr())
library(readr)
library(dplyr)
#Use the library rjson to read the JSON file.
library(rjson)


#################
# Funciones
#################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}












thr<-0.8

file92<- paste("Kink_NT_Map92_",thr,".csv", sep = "")
file93<- paste("Kink_NT_Map93_",thr,".csv", sep = "")
file94<- paste("Kink_NT_Map94_",thr,".csv", sep = "")
file95<- paste("Kink_NT_Map95_",thr,".csv", sep = "")

file100<- paste("Kink_CT_Map100_",thr,".csv", sep = "")
file101<- paste("Kink_CT_Map101_",thr,".csv", sep = "")
file102<- paste("Kink_CT_Map102_",thr,".csv", sep = "")
file103<- paste("Kink_CT_Map103_",thr,".csv", sep = "")


angles92 <- read_delim(file92, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles92<-angles92[,1:2]
angles93 <- read_delim(file93, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles93<-angles93[,1:2]
angles94 <- read_delim(file94, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles94<-angles94[,1:2]
angles95 <- read_delim(file95, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles95<-angles95[,1:2]




 
angles100 <- read_delim(file100, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles100<-angles100[,1:2]
angles101 <- read_delim(file101, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles101<-angles101[,1:2]
angles102 <- read_delim(file102, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles102<-angles102[,1:2]
angles103 <- read_delim(file103, delim = ";", escape_double = FALSE, trim_ws = TRUE)
angles103<-angles103[,1:2]










############### 92


map92 <- ggplot(angles92, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 92")+
  geom_bin_2d(bins = 72,aes(fill = ..ndensity..)) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(map92)




linesmap92 <- ggplot(angles92, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 92")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap92)




pointsmap92 <- ggplot(angles92, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 92")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap92)





############  93



map93 <- ggplot(angles93, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 93")+
  geom_bin_2d(bins = 72,aes(fill = ..ndensity..)) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(map93)




linesmap93 <- ggplot(angles93, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 93")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap93)


pointsmap93 <- ggplot(angles93, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 93")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap93)





############  94



map94 <- ggplot(angles94, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 94")+
  geom_bin_2d(bins = 72,show.legend = FALSE) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(map94)




linesmap94 <- ggplot(angles94, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 94")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap94)


pointsmap94 <- ggplot(angles94, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 94")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap94)




############  95



map95 <- ggplot(angles95, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 95")+
  geom_bin_2d(bins = 72,show.legend = FALSE) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(map95)




linesmap95 <- ggplot(angles95, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 95")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap95)


pointsmap95 <- ggplot(angles95, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 95")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap95)

















################# 100


map100 <- ggplot(angles100, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 100")+
  geom_bin_2d(bins = 72,show.legend = FALSE) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


  #theme(legend.title = element_text(size = 15),legend.position="top")
print(map100)




linesmap100 <- ggplot(angles100, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 100")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap100)




pointsmap100 <- ggplot(angles100, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 100")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap100)









############  101



map101 <- ggplot(angles101, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 101")+
  geom_bin_2d(bins = 72,show.legend = FALSE) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(map101)




linesmap101 <- ggplot(angles101, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 101")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap101)


pointsmap101 <- ggplot(angles101, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 101")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap101)





############  102



map102 <- ggplot(angles102, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 102")+
  geom_bin_2d(bins = 72,show.legend = FALSE) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(map101)




linesmap102 <- ggplot(angles102, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 102")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap102)


pointsmap102 <- ggplot(angles102, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 102")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap102)






############  103



map103 <- ggplot(angles103, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 103")+
  geom_bin_2d(bins = 72,show.legend = FALSE) +
  xlim(c(-185, 185)) + ylim(c(-185, 185)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(map103)




linesmap103 <- ggplot(angles103, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 103")+
  geom_density2d(bins = 20)+
  xlim(c(-181, 181)) + ylim(c(-181, 181)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(linesmap103)


pointsmap103 <- ggplot(angles103, aes(x = PHI, y = PSI)) +
  labs(title="Mapa de ramachandran 103")+
  geom_jitter()+
  xlim(c(-180, 180)) + ylim(c(-180, 180)) +
  #scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  #scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  geom_hline(yintercept=180)+ geom_hline(yintercept=-180)+
  geom_vline(xintercept=180)+ geom_vline(xintercept=-180)+
  theme(legend.position = "none")+
  theme_bw()


#theme(legend.title = element_text(size = 15),legend.position="top")
print(pointsmap103)





################## multiplots ##############3

listN <- list(map92,map93,map94,map95,linesmap92,linesmap93,linesmap94,linesmap95,pointsmap92,pointsmap93,pointsmap94,pointsmap95)

matN<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=3, byrow=TRUE)
multiplot(plotlist=listN, layout=matN)




listC <- list(map100,map101,map102,map103,linesmap100,linesmap101,linesmap102,linesmap103,pointsmap100,pointsmap101,pointsmap102,pointsmap103)

matC<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=3, byrow=TRUE)
#multiplot(plotlist=listC, layout=matC)





############################################################################
############################################################################
############################################################################


#files direction: /home/ivan/disk/PhD/loops_H3/kink_sequence/protes



map92_1<-read.csv("Mapa_92_1.txt",sep = ";")




map92_1plot <- ggplot(map92_1,aes(x=PHI, y=PSI, fill=VALUE)) +
  labs(title="Mapa de 192")+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+
  scale_x_continuous(breaks = seq(-180, 180, by = 60))+
  scale_y_continuous(breaks = seq(-180, 180, by = 60))+
  theme_bw()

print(map92_1plot)










