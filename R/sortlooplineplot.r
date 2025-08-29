#Ivan Martin Hernandez

#Scrip para la representacion de la figura de la mejora de Korpe respecto al paper original.
library(ggplot2) 

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



pintame =function(matriz){
  
  tt<- as.data.frame(t(matriz))
  di<-dim(matriz)[2]
  time<-1:di
  tt <- cbind(tt, time)
  
  tt <- melt(tt ,  id.vars = 'time', variable.name = 'series')
  print(ggplot(tt, aes(time,value)) + geom_line(aes(colour = series)))
}




#################
# carga de datos
#################


# carga de los nombres de los pdb
names<- c("1mlb", "1x9q", "2d7t", "2e27", "3g5y", "3hc4", "1jpt","3e8u", "3m8o", "1mfa", "1mqk", "1nlb", "2adf", "2fbj", "2w60", "3gnm", "3hnt", "3v0w", "1dlf", "2xwt", "2ypv", "3ifl", "3liz", "3mxw", "3oz9", "3umt", "4h0h", "4h20", "1oaq", "2v17", "3t65", "4hpy", "1jfq", "2r8s", "2vxv", "3eo9", "3i9g", "3p0y", "3giz", "1fns", "1gig", "1seq", "3go1", "3mlr", "3lmj", "4f57", "4nzu", "2fb4", "3nps")
size<- c(9, 9, 9, 9, 9, 9, 10,10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 16, 16, 16, 16, 18, 18, 18, 18, 19, 19)
Dataestandar<-data.frame(Real_Loop=names,  size=size)

colores<- c("#f75a00", "#b4f733", "#0043f7", "#43e633", "#7e261f", "#f0c824", "#e400f7", "#30bec4", "#dc1c1c", "#6c6464")
#colores<-c("#f75a00", "#b4f733", "#0043f7", "#43e633", "#e400f7", "#30bec4", "#f7d200", "#2f63b3", "#dc1c1c", "#6c6464")



### NATIVE-CRYSTAL
rosetta1_table <- read.table("./Kink_improve_50K/A_rosetta1_native_50K_k090160.txt", sep = "", header = FALSE)
rosetta2_table <- read.table("./Kink_improve_50K/A_rosetta2_native_50K_k090160.txt", sep = "", header = FALSE)
rosetta5_table <- read.table("./Kink_improve_50K/A_rosetta5_native_50K_k090160.txt", sep = "", header = FALSE)
rosetta10_table <- read.table("./Kink_improve_50K/A_rosetta10_native_50K_k090160.txt", sep = "", header = FALSE)
#rosetta10_table <- read.table("./Kink_improve_50K/A_rosetta10_modeling_50K_k090160.txt", sep = "", header = FALSE)
#rosetta10_table <- read.table("./Kink_improve_50K/rosetta10_mod_modeling_50K_k090160.txt", sep = "", header = FALSE)


DMtop1<-rosetta1_table[,c(1,4,5,9,12)]
colnames(DMtop1)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop1["Loop"],  value=DMtop1["RMSD_BBO"])

#maxvalue<- max(DMtop1["RMSD_BBO"])
maxvalue<- 10.4

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop1 <- merge(data,Dataestandar,by="Real_Loop")

a <- ggplot(FinalTop1, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-native Top 1", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))


DMtop2<-rosetta2_table[,c(1,4,5,9,12)]
colnames(DMtop2)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop2["Loop"],  value=DMtop2["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop2 <- merge(data,Dataestandar,by="Real_Loop")

b <- ggplot(FinalTop2, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-native Top 2", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))


DMtop5<-rosetta5_table[,c(1,4,5,9,12)]
colnames(DMtop5)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop5["Loop"],  value=DMtop5["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop5 <- merge(data,Dataestandar,by="Real_Loop")

c <- ggplot(FinalTop5, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-native Top 5", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))



DMtop10<-rosetta10_table[,c(1,4,5,9,12)]
colnames(DMtop10)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop10["Loop"],  value=DMtop10["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop10 <- merge(data,Dataestandar,by="Real_Loop")

d <- ggplot(FinalTop10, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-native Top 10", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =13), axis.title = element_text(size = 12), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.text.y = element_text(size = 10),legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))


zz <- ggplot(FinalTop10, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-native Top 10", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 15), legend.title = element_text(size = 16),legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(0.8, "cm"))

plot(zz)


### MODELING-CRYSTAL
rosetta1_table <- read.table("./Kink_improve_50K/A_rosetta1_modeling_50K_k090160.txt", sep = "", header = FALSE)
rosetta2_table <- read.table("./Kink_improve_50K/A_rosetta2_modeling_50K_k090160.txt", sep = "", header = FALSE)
rosetta5_table <- read.table("./Kink_improve_50K/A_rosetta5_modeling_50K_k090160.txt", sep = "", header = FALSE)
rosetta10_table <- read.table("./Kink_improve_50K/A_rosetta10_modeling_50K_k090160.txt", sep = "", header = FALSE)

DMtop1<-rosetta1_table[,c(1,4,5,9,12)]
colnames(DMtop1)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop1["Loop"],  value=DMtop1["RMSD_BBO"])

#maxvalue<- max(DMtop1["RMSD_BBO"])
maxvalue<- 10.4

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop1 <- merge(data,Dataestandar,by="Real_Loop")

e <- ggplot(FinalTop1, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-modeling Top 1", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))


DMtop2<-rosetta2_table[,c(1,4,5,9,12)]
colnames(DMtop2)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop2["Loop"],  value=DMtop2["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop2 <- merge(data,Dataestandar,by="Real_Loop")

f <- ggplot(FinalTop2, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-modeling Top 2", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))


DMtop5<-rosetta5_table[,c(1,4,5,9,12)]
colnames(DMtop5)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop5["Loop"],  value=DMtop5["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop5 <- merge(data,Dataestandar,by="Real_Loop")

g <- ggplot(FinalTop5, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-modeling Top 5", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))



DMtop10<-rosetta10_table[,c(1,4,5,9,12)]
colnames(DMtop10)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop10["Loop"],  value=DMtop10["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop10 <- merge(data,Dataestandar,by="Real_Loop")

h <- ggplot(FinalTop10, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Crystal-modeling Top 10", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =13), axis.title = element_text(size = 12), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.text.y = element_text(size = 10),legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))

print(h)


### MODELING-MODEL
rosetta1_table <- read.table("./Kink_improve_50K/rosetta1_mod_modeling_50K_k090160.txt", sep = "", header = FALSE)
rosetta2_table <- read.table("./Kink_improve_50K/rosetta2_mod_modeling_50K_k090160.txt", sep = "", header = FALSE)
rosetta5_table <- read.table("./Kink_improve_50K/rosetta5_mod_modeling_50K_k090160.txt", sep = "", header = FALSE)
rosetta10_table <- read.table("./Kink_improve_50K/rosetta10_mod_modeling_50K_k090160.txt", sep = "", header = FALSE)


DMtop1<-rosetta1_table[,c(1,4,5,9,12)]
colnames(DMtop1)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop1["Loop"],  value=DMtop1["RMSD_BBO"])

#maxvalue<- max(DMtop1["RMSD_BBO"])
maxvalue<- 10.4

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop1 <- merge(data,Dataestandar,by="Real_Loop")

z <- ggplot(FinalTop1, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Model-modeling Top 1", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))


DMtop2<-rosetta2_table[,c(1,4,5,9,12)]
colnames(DMtop2)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop2["Loop"],  value=DMtop2["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop2 <- merge(data,Dataestandar,by="Real_Loop")

j <- ggplot(FinalTop2, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Model-modeling Top 2", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))


DMtop5<-rosetta5_table[,c(1,4,5,9,12)]
colnames(DMtop5)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop5["Loop"],  value=DMtop5["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop5 <- merge(data,Dataestandar,by="Real_Loop")

k <- ggplot(FinalTop5, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Model-modeling Top 5", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =11), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))



DMtop10<-rosetta10_table[,c(1,4,5,9,12)]
colnames(DMtop10)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
data=data.frame(name=DMtop10["Loop"],  value=DMtop10["RMSD_BBO"])

for (i in 1:49){
  pdb= substr(data[i,1], 1, 4)
  data[i,"Real_Loop"]<-pdb
}

FinalTop10 <- merge(data,Dataestandar,by="Real_Loop")


l <- ggplot(FinalTop10, aes(x=Real_Loop, y=RMSD_BBO)) + 
  geom_bar(stat = "identity", aes(fill=as.factor(size)))+
  xlim(names) +
  scale_fill_manual(values = colores ) +
  labs(title="Model-modeling Top 10", x=element_blank(), y="RMSD", fill="Loop size") +
  ylim(c(NA, maxvalue))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5,size =13), axis.title = element_text(size = 12), axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.text.y = element_text(size = 10),legend.title = element_text(size = 9),legend.text = element_text(size = 8, face = "bold"),legend.key.size = unit(0.4, "cm"))

plot(l)


#plot the grap

tt<-list(a,b,c,d,e,f,g,h,z,j,k,l)

pp<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=4, byrow=FALSE)

multiplot(plotlist=tt, layout=pp)

final<-list (d,h,l)

pp<-matrix(c(1,2,3), nrow=3, byrow=FALSE)

multiplot(plotlist=final, layout=pp)
