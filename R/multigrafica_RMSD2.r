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
#names <- sort(names)

multilistor <-list()
multilist <-list()
multilist2 <-list()

for (j in 1:length(names)){
  #totalname<- paste(names[j],"_5_modeling_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_modeling_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  
  totalname<- paste(names[j],"_4_native_50K_k090160_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_modeling_50K_k090160_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_5_modeling_50K_k090160_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_native_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  
  #totalname<- paste(names[j],"_4_native_newminimization_rosettaS.txt", sep = "")

  if(file.exists(totalname)){
    rosetta_table <- read.table(totalname, sep = "", header = FALSE)
    DM1<-rosetta_table[,c(1,4,5,9,12)]
    colnames(DM1)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
    
    DM2<- DM1[sort.list(DM1[,4]), ]
    DM3<-subset(DM2[1:800,])
    
    Average<-mean(DM1[,4])
    sd<-sd(DM1[,4])
    ylimit<-Average + (2*sd)
    
    Average2<-mean(DM3[,4])
    sd2<-sd(DM3[,4])
    ylimit2<-Average2 + (2.5*sd2)
    
    for (i in 1:1000){
      if (i<=10){DM2[i,"grupo"]<-"top10"}
      else{DM2[i,"grupo"]<-"notop"}
    }
    
    a <- ggplot(DM2, aes(x=RMSD_BBO, y=ERosetta, color=grupo) ) + 
      geom_point( aes(color=grupo, size=grupo)) +
      scale_size_manual(values=c(0,1)) +
      scale_color_manual(values=c('#999999','#E69F00')) +
      xlim(c(-0, 10)) +
      labs(title=names[j],x=element_blank(), y=element_blank()) +
      theme_classic()+
      theme(legend.position = "none")+
      theme(plot.title = element_text(hjust=0.5,size = 9))
    
    multilistor[j]<-list(a)
    
    a<-a+ylim(c(NA, ylimit))
    multilist[j]<-list(a)
    
    a<-a+ylim(c(NA, ylimit2))
    multilist2[j]<-list(a)
    
  }
  else{
    #print("No existe")
    b <- ggplot(DM2, aes(x=RMSD_BBO, y=ERosetta) ) + 
      geom_point(size=-1) +
      xlim(c(-0, 10)) +
      labs(title=names[j],x=element_blank(), y=element_blank()) +
      theme_classic()+
      theme(plot.title = element_text(hjust=0.5,size = 9))
    
    multilistor[j]<-list(b)
    multilist[j]<-list(b)
    multilist2[j]<-list(b)
    }
}

pp<-matrix(c(1,2,3,4,5,6,7,8,9,"",11,12), nrow=4, byrow=FALSE)
pp<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49), nrow=7, byrow=TRUE)
#pp<-matrix(c(1,2,3,4,5,6,"",7,8,9,10,11,12,13,14,15,16,17,18,"",19,20,21,22,23,24,25,26,27,28,29,30,31,32,"","","","","","",33,34,35,36,37,38,"","",39,"",40,41,42,43,"","","","","","",44,45,46,47,"","","","",48,49), nrow=7, byrow=TRUE)

#multiplot(plotlist=multilist, cols = 7)
multiplot(plotlist=multilist2, layout=pp)
multiplot(plotlist=multilistor, cols = 7)


multilistkorp <-list()

for (j in 1:length(names)){
  #totalname<- paste(names[j],"_5_modeling_50K_k090160_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_modeling_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_native_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_modeling_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  
  
  if(file.exists(totalname)){
    rosetta_table <- read.table(totalname, sep = "", header = FALSE)
    DM1<-rosetta_table[,c(1,4,5,9,12)]
    colnames(DM1)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
    
    DM2<- DM1[sort.list(DM1[,4]), ]
    DM3<-subset(DM2[1:800,])
    
    Average<-mean(DM1[,4])
    sd<-sd(DM1[,4])
    ylimit<-Average + (2*sd)
    
    Average2<-mean(DM3[,4])
    sd2<-sd(DM3[,4])
    ylimit2<-Average2 + (2.5*sd2)
    
    for (i in 1:1000){
      if (i<=10){DM2[i,"grupo"]<-"top10"}
      else{DM2[i,"grupo"]<-"notop"}
    }
    
    a <- ggplot(DM2, aes(x=RMSD_RCD, y=EKORP, color=grupo) ) + 
      geom_point( aes(color=grupo, size=grupo)) +
      scale_size_manual(values=c(0,1)) +
      scale_color_manual(values=c('#999999','#E69F00')) +
      xlim(c(-0, 10)) +
      labs(title=names[j],x=element_blank(), y=element_blank()) +
      theme_classic()+
      theme(legend.position = "none")+
      theme(plot.title = element_text(hjust=0.5,size = 9))
    
    multilistkorp[j]<-list(a)
    
    
  }
  else{
    #print("No existe")
    b <- ggplot(DM1, aes(x=RMSD_RCD, y=EKORP) ) + 
      geom_point(size=-1) +
      xlim(c(-0, 10)) +
      labs(title=names[j],x=element_blank(), y=element_blank()) +
      theme_classic()+
      theme(plot.title = element_text(hjust=0.5,size = 9))
    
    multilistkorp[j]<-list(b)
  }
}

multiplot(plotlist=multilistkorp, cols = 7)







rosetta_table <- read.table("3g5y_4_native_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "", header = FALSE)


DM1<-rosetta_table[,c(1,4,5,9,12)]
colnames(DM1)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")

DM2<- DM1[sort.list(DM1[,4]), ]
DM3<-subset(DM2[1:800,])

Average<-mean(DM1[,4])
sd<-sd(DM1[,4])
ylimit<-Average + (2*sd)

Average2<-mean(DM3[,4])
sd2<-sd(DM3[,4])
ylimit2<-Average2 + (2.5*sd2)


for (i in 1:1000){
  if (i<=10){DM1[i,"grupo"]<-"top10"}
  else{DM1[i,"grupo"]<-"notop"}
}

a <- ggplot(DM1, aes(x=RMSD_BBO, y=ERosetta, color=grupo) ) + 
  #geom_point(size=0) +
  geom_point( aes(color=grupo, size=grupo)) +
  scale_size_manual(values=c(0,1)) +
  scale_color_manual(values=c('#999999','#E69F00')) +
  xlim(c(-0, 10)) +
  labs(title="3umt",x=element_blank(), y=element_blank()) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust=0.5,size = 9))
  #ylim(c(NA, ylimit))
  
plot(a)

b <- ggplot(DM1, aes(x=RMSD_BBO, y=ERosetta) ) + 
  geom_point(size=-1) +
  xlim(c(-0, 10)) +
  theme_bw()


plot(b)

tt<-list(a,b,a,a,a,a,a,a,a)

multiplot(plotlist=tt, cols = 3)























# carga de los nombres de los pdb
names<- c("1dlf", "2ypv", "4nzu")
#names <- sort(names)

multilistor <-list()
multilist <-list()
multilist2 <-list()



for (j in 1:length(names)){
  #totalname<- paste(names[j],"_5_modeling_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_modeling_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_native_rcd_H3P_n5Hk1k99t_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_4_native_50K_k090160_rosettaS.txt", sep = "")
  totalname<- paste(names[j],"_4_modeling_50K_k090160_rosettaS.txt", sep = "")
  #totalname<- paste(names[j],"_5_modeling_50K_k090160_rosettaS.txt", sep = "")
  
  
  
  if(file.exists(totalname)){
    rosetta_table <- read.table(totalname, sep = "", header = FALSE)
    DM1<-rosetta_table[,c(1,4,5,9,12)]
    colnames(DM1)<-c("Loop", "RMSD_RCD","EKORP","ERosetta", "RMSD_BBO")
    
    DM2<- DM1[sort.list(DM1[,4]), ]
    DM3<-subset(DM2[1:800,])
    
    Average<-mean(DM1[,4])
    sd<-sd(DM1[,4])
    ylimit<-Average + (2*sd)
    
    Average2<-mean(DM3[,4])
    sd2<-sd(DM3[,4])
    ylimit2<-Average2 + (2.5*sd2)
    
    for (i in 1:1000){
      if (i<=10){DM2[i,"grupo"]<-"top10"}
      else{DM2[i,"grupo"]<-"notop"}
    }
    
    a <- ggplot(DM2, aes(x=RMSD_BBO, y=ERosetta, color=grupo) ) + 
      geom_point( aes(color=grupo, size=grupo)) +
      scale_size_manual(values=c(0,1)) +
      scale_color_manual(values=c('#999999','#E69F00')) +
      xlim(c(-0, 10)) +
      labs(title=names[j],x="RMSD", y="Rosetta energy" ) +
      theme_classic()+
      theme(legend.position = "none")+
      theme(plot.title = element_text(hjust=0.5,size = 14), axis.title = element_text(size = 12),axis.text.y = element_text(size = 10))
    
    multilistor[j]<-list(a)
    
    a<-a+ylim(c(NA, ylimit))
    multilist[j]<-list(a)
    
    a<-a+ylim(c(NA, ylimit2))
    multilist2[j]<-list(a)
    
  }
  else{
    #print("No existe")
    b <- ggplot(DM2, aes(x=RMSD_BBO, y=ERosetta) ) + 
      geom_point(size=-1) +
      xlim(c(-0, 10)) +
      labs(title=names[j],x=element_blank(), y=element_blank()) +
      theme_classic()+
      theme(plot.title = element_text(hjust=0.5,size = 14), axis.title = element_text(size = 12))
    
    multilistor[j]<-list(b)
    multilist[j]<-list(b)
    multilist2[j]<-list(b)
  }
}

pp<-matrix(c(1,2,3), nrow=1, byrow=FALSE)

#multiplot(plotlist=multilist, cols = 7)
multiplot(plotlist=multilist2, layout=pp)
multiplot(plotlist=multilistor, cols = 7)





