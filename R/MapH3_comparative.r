#20/11/2022   Ivan Martin Hernandez
#Script that generate the figures of RCD H3 loops with/out maps


#Libraries
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readr)
theme_set(theme_pubr())
library(reshape2)
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






#################
# Main
#################


############# Modeling ########

fileNKrep1<-"Modeling/Rep1/AVg_modeling_50K_2022_nK.txt"
fileNKrep2<-"Modeling/Rep2/AVg_modeling_50K_2022_nK.txt"
fileNKrep3<-"Modeling/Rep1/AVg_modeling_50K_2022_nK.txt"
fileKrep1<-"Modeling/Rep1/AVg_modeling_50K_2022S_90_070_3.txt"
fileKrep2<-"Modeling/Rep2/AVg_modeling_50K_2022S_90_070_3.txt"
fileKrep3<-"Modeling/Rep3/AVg_modeling_50K_2022S_90_070_3.txt"


tableNKrep1 <- read_delim(fileNKrep1, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableNKrep2 <- read_delim(fileNKrep2, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableNKrep3 <- read_delim(fileNKrep3, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)

tableKrep1 <- read_delim(fileKrep1, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableKrep2 <- read_delim(fileKrep2, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableKrep3 <- read_delim(fileKrep3, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)



NKtable <- bind_cols(tableNKrep1[,37],tableNKrep2[,37],tableNKrep3[,37])
colnames(NKtable) <- c("Rep1","Rep2","Rep3")
NKtable <- as.data.frame(NKtable)
NKtable <- melt(NKtable)


Ktable <- bind_cols(tableKrep1[,35],tableKrep2[,36],tableKrep3[,36])
colnames(Ktable) <- c("Rep1","Rep2","Rep3")
Ktable <- as.data.frame(Ktable)
Ktable <- melt(Ktable)

Ttable <- bind_cols(tableNKrep1[,37],tableNKrep2[,37],tableNKrep3[,37],tableKrep1[,35],tableKrep2[,36],tableKrep3[,36])
colnames(Ttable) <- c("Sin Mapa Rep1","Sin Mapa Rep2","Sin Mapa Rep3","Con Mapa Rep1","Con Mapa Rep2","Con Mapa Rep3")
Ttable <- as.data.frame(Ttable)
Ttable <- melt(Ttable)

Ttable$modo = c("Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa")

############# plots ########


NKboxplotM <- ggplot(data = NKtable, aes(x=variable, y=value)) +
  labs(title="Replicas Modeling sin Mapas H3")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(NKboxplotM)

NboxplotM <- ggplot(data = Ktable, aes(x=variable, y=value)) +
  labs(title="Replicas Modeling con Mapas H3")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(NboxplotM)

TboxplotM <- ggplot(data = Ttable, aes(x=modo, y=value)) +
  labs(title="Replicas Modeling")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  guides(fill=guide_legend(title="Replicas"))+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  #theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  theme( plot.title = element_text(hjust = 0.5))

print(TboxplotM)

TyboxplotM <- ggplot(data = Ttable, aes(x=modo, y=value)) +
  labs(title="Resultados Modeling")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=modo))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  #theme( plot.title = element_text(hjust = 0.5))

print(TyboxplotM)




############# Modeling HQ ########

fileNKrep1<-"Modeling_HQ/Rep1/AVg_modeling_50K_2022_nK.txt"
fileNKrep2<-"Modeling_HQ/Rep2/AVg_modeling_50K_2022_nK.txt"
fileNKrep3<-"Modeling_HQ/Rep1/AVg_modeling_50K_2022_nK.txt"
fileKrep1<-"Modeling_HQ/Rep1/AVg_modeling_50K_2022S_90_070_3.txt"
fileKrep2<-"Modeling_HQ/Rep2/AVg_modeling_50K_2022S_90_070_3.txt"
fileKrep3<-"Modeling_HQ/Rep3/AVg_modeling_50K_2022S_90_070_3.txt"


tableNKrep1 <- read_delim(fileNKrep1, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableNKrep2 <- read_delim(fileNKrep2, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableNKrep3 <- read_delim(fileNKrep3, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)

tableKrep1 <- read_delim(fileKrep1, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableKrep2 <- read_delim(fileKrep2, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableKrep3 <- read_delim(fileKrep3, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)



NKtable <- bind_cols(tableNKrep1[,36],tableNKrep2[,36],tableNKrep3[,36])
colnames(NKtable) <- c("Rep1","Rep2","Rep3")
NKtable <- as.data.frame(NKtable)
NKtable <- melt(NKtable)


Ktable <- bind_cols(tableKrep1[,36],tableKrep2[,36],tableKrep3[,36])
colnames(Ktable) <- c("Rep1","Rep2","Rep3")
Ktable <- as.data.frame(Ktable)
Ktable <- melt(Ktable)

Ttable <- bind_cols(tableNKrep1[,36],tableNKrep2[,36],tableNKrep3[,36],tableKrep1[,36],tableKrep2[,36],tableKrep3[,36])
colnames(Ttable) <- c("Sin Mapa Rep1","Sin Mapa Rep2","Sin Mapa Rep3","Con Mapa Rep1","Con Mapa Rep2","Con Mapa Rep3")
Ttable <- as.data.frame(Ttable)
Ttable <- melt(Ttable)

Ttable$modo = c("Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa")

############# plots ########


NKboxplotHQ <- ggplot(data = NKtable, aes(x=variable, y=value)) +
  labs(title="Replicas Modeling HQ sin Mapas H3")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(NKboxplotHQ)

NboxplotHQ <- ggplot(data = Ktable, aes(x=variable, y=value)) +
  labs(title="Replicas Modeling HQ con Mapas H3")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(NboxplotHQ)

TboxplotHQ <- ggplot(data = Ttable, aes(x=modo, y=value)) +
  labs(title="Replicas Modeling HQ")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  guides(fill=guide_legend(title="Replicas"))+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  #theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  theme( plot.title = element_text(hjust = 0.5))

print(TboxplotHQ)

TyboxplotHQ <- ggplot(data = Ttable, aes(x=modo, y=value)) +
  labs(title="Resultados Modeling HQ")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=modo))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
#theme( plot.title = element_text(hjust = 0.5))

print(TyboxplotHQ)






############# Native ########

fileNKrep1<-"Nativo/Rep1/AVg_native_50K_2022_nK.txt"
fileNKrep2<-"Nativo/Rep2/AVg_native_50K_2022_nK.txt"
fileNKrep3<-"Nativo/Rep1/AVg_native_50K_2022_nK.txt"
fileKrep1<-"Nativo/Rep1/AVg_native_50K_2022S_90_070_3.txt"
fileKrep2<-"Nativo/Rep2/AVg_native_50K_2022S_90_070_3.txt"
fileKrep3<-"Nativo/Rep3/AVg_native_50K_2022S_90_070_3.txt"


tableNKrep1 <- read_delim(fileNKrep1, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableNKrep2 <- read_delim(fileNKrep2, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableNKrep3 <- read_delim(fileNKrep3, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)

tableKrep1 <- read_delim(fileKrep1, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableKrep2 <- read_delim(fileKrep2, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)
tableKrep3 <- read_delim(fileKrep3, delim = " ",col_names = FALSE,escape_double = FALSE, trim_ws = FALSE)



NKtable <- bind_cols(tableNKrep1[,38],tableNKrep2[,36],tableNKrep3[,38])
colnames(NKtable) <- c("Rep1","Rep2","Rep3")
NKtable <- as.data.frame(NKtable)
NKtable <- melt(NKtable)


Ktable <- bind_cols(tableKrep1[,36],tableKrep2[,36],tableKrep3[,36])
colnames(Ktable) <- c("Rep1","Rep2","Rep3")
Ktable <- as.data.frame(Ktable)
Ktable <- melt(Ktable)

Ttable <- bind_cols(tableNKrep1[,38],tableNKrep2[,36],tableNKrep3[,38],tableKrep1[,36],tableKrep2[,36],tableKrep3[,36])
colnames(Ttable) <- c("Sin Mapa Rep1","Sin Mapa Rep2","Sin Mapa Rep3","Con Mapa Rep1","Con Mapa Rep2","Con Mapa Rep3")
Ttable <- as.data.frame(Ttable)
Ttable <- melt(Ttable)

Ttable$modo = c("Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Sin Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa","Con Mapa")

############# plots ########


NKboxplotN <- ggplot(data = NKtable, aes(x=variable, y=value)) +
  labs(title="Replicas Nativo sin Mapas H3")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(NKboxplotN)

NboxplotN <- ggplot(data = Ktable, aes(x=variable, y=value)) +
  labs(title="Replicas Nativo con Mapas H3")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(NboxplotN)

TboxplotN <- ggplot(data = Ttable, aes(x=modo, y=value)) +
  labs(title="Replicas Nativo")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  guides(fill=guide_legend(title="Replicas"))+
  geom_boxplot(aes(fill=variable))+
  ylim(c(0,3.05)) +
  theme_bw()+
  #theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  theme( plot.title = element_text(hjust = 0.5))

print(TboxplotN)

TyboxplotN <- ggplot(data = Ttable, aes(x=modo, y=value)) +
  labs(title="Resultados Nativo")+
  ylab(label="RMSD (Å)")+ xlab(label=NULL)+
  geom_boxplot(aes(fill=modo))+
  ylim(c(0,3.05)) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
#theme( plot.title = element_text(hjust = 0.5))

print(TyboxplotN)





#################
# Multiplots
#################

list <- list(NKboxplotM,NKboxplotHQ,NKboxplotN,NboxplotM,NboxplotHQ,NboxplotN,TboxplotM,TboxplotHQ,TboxplotN,TyboxplotM,TyboxplotHQ,TyboxplotN)
matN<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=4, byrow=TRUE)

multiplot(plotlist=list, layout=matN)










