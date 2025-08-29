

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


df1 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.1S.txt", header=FALSE)

subdf1<-df1$V10

a <- ggplot(df1, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.1    size=",length(df1$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(a)



df2 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.2S.txt", header=FALSE)
b <- ggplot(df2, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.2    size=",length(df2$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(b)


df3 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.3S.txt", header=FALSE)
c <- ggplot(df3, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.3    size=",length(df3$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(c)

df4 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.4S.txt", header=FALSE)
d <- ggplot(df4, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.4    size=",length(df4$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(d)

df5 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.5S.txt", header=FALSE)
e <- ggplot(df5, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.5    size=",length(df5$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(e)

df6 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.6S.txt", header=FALSE)
f <- ggplot(df6, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.6    size=",length(df6$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(f)

df7 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.7S.txt", header=FALSE)
g <- ggplot(df7, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.7    size=",length(df7$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(g)

df8 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.8S.txt", header=FALSE)
h <- ggplot(df8, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.8    size=",length(df8$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(h)

df9 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC0.9S.txt", header=FALSE)
i <- ggplot(df9, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 0.9    size=",length(df9$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(i)


df10 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC1S.txt", header=FALSE)
j <- ggplot(df10, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 1.0    size=",length(df10$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(j)

df11 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_lazyC1.2S.txt", header=FALSE)
k <- ggplot(df11, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Lazy 1.2    size=",length(df11$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(k)

df12 <- read.table("3umt_4_native_ori_50K_3umt_2022S_iv_95_070_3_rosetta_RClustS.txt", header=FALSE)
l <- ggplot(df12, aes(x=V10)) + 
  geom_density()+
  labs(title= paste("Rclust    size=",length(df12$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(l)



tt<-list(a,b,c,d,e,f,g,h,i,j,k,l)

pp<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=4, byrow=FALSE)

multiplot(plotlist=tt, layout=pp)


#before Rosetta


a <- ggplot(df1, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.1    size=",length(df1$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(a)

b <- ggplot(df2, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.2    size=",length(df2$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(b)


c <- ggplot(df3, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.3    size=",length(df3$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(c)

d <- ggplot(df4, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.4    size=",length(df4$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(d)

e <- ggplot(df5, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.5    size=",length(df5$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(e)

f <- ggplot(df6, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.6    size=",length(df6$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(f)

g <- ggplot(df7, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.7    size=",length(df7$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(g)

h <- ggplot(df8, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.8    size=",length(df8$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(h)

i <- ggplot(df9, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 0.9    size=",length(df9$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(i)


j <- ggplot(df10, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 1.0    size=",length(df10$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(j)

k <- ggplot(df11, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Lazy 1.2    size=",length(df11$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(k)

l <- ggplot(df12, aes(x=V6)) + 
  geom_density()+
  labs(title= paste("Rclust    size=",length(df12$V10) , sep = " "), x="RMSD") +
  xlim(c(0, 8))
plot(l)



tt<-list(a,b,c,d,e,f,g,h,i,j,k,l)

pp<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=4, byrow=FALSE)

multiplot(plotlist=tt, layout=pp)


