## ----load libraries, echo=FALSE-----------------------------------------------
library(ecotraj)
library(tidyr) ## For data manipulation
library(ggplot2) ## For plotting
library(hrbrthemes) ## Additional themes for ggplot2
library(scales) ## Scale functions for visualisation
library(viridis) ## Viridis color map

## ----load furseals, echo=T----------------------------------------------------
data("furseals")

## ---- echo=T------------------------------------------------------------------
Net_changes<-trajectoryLengths2D(furseals[,c("d13C","d15N")],
                                 furseals$ID_SITA,
                                 furseals$Time, relativeToInitial=TRUE) 
head(Net_changes)

## ---- echo=T------------------------------------------------------------------
Segment_lengths<-trajectoryLengths2D(furseals[,c("d13C","d15N")],
                                     furseals$ID_SITA,
                                     furseals$Time, relativeToInitial=FALSE) 
head(Segment_lengths)

## ---- echo=T------------------------------------------------------------------
Angles<-trajectoryAngles2D(furseals[,c("d13C","d15N")],
                           furseals$ID_SITA,
                           furseals$Time, betweenSegments=FALSE)
head(Angles)

## ---- echo=T------------------------------------------------------------------
D <- dist(furseals[,c("d13C","d15N")])
Ds<-trajectoryDistances(D, furseals$ID_SITA, surveys = NULL, distance.type = "DSPD",
                        symmetrization = "mean", add = TRUE, verbose = FALSE)

## ---- echo=TRUE, fig.height=4, fig.width=6------------------------------------
colstd<-c("black","yellow","green","blue","grey","red")
pt<-c(16,16,16,16)
hsxy <- hclust(Ds, "ward.D2")
plot(hsxy,hang = -1, main="distance Fur Seals", cex=.6)
Hst=2 # Cutting height
x<-rect.hclust(hsxy, h=Hst,
               border = colstd)

## -----------------------------------------------------------------------------
groups <- cutree(hsxy, h=Hst)
furseals$cluster <- as.factor(groups)

## ---- echo=T------------------------------------------------------------------
furseals$sp_gender<-paste(furseals$Sexe, furseals$Species, sep=" ")

## ---- echo=TRUE, fig.height=6, fig.width=6------------------------------------
ggplot(data=furseals,aes(x=d13C,y=d15N,color=cluster,shape=Place))+
  geom_point()+
  geom_path(aes(x=d13C,y=d15N,group=ID_SITA,color=cluster),
            arrow = arrow(length = unit(0.10, "cm")))+
  xlab(expression(delta^13*"C"))+
  ylab(expression(delta^15*"N"))+
  facet_wrap(~sp_gender) +
  theme_classic()

## -----------------------------------------------------------------------------
NC<-Net_changes[,-30]
NC$cluster<-furseals$cluster[1:47]
NC$ID<-as.numeric(rownames(NC))
colnames(NC)<-c(2:30,"cluster","ID")

## -----------------------------------------------------------------------------
NCline <- tidyr::pivot_longer(NC, 1:29, 
                              names_to ="Time_from_present", 
                              values_to="Net_changes", 
                              names_transform = function(x) {as.numeric(x)-1}) 
colnames(NCline)[1:2]<-c("Clusters", "ID")
NCline <- NCline[order(NCline$Time_from_present, decreasing=FALSE),]
NCline <- as.data.frame(NCline)
NCline$sp_gender<-c(furseals$sp_gender[1:47])

## ---- fig.height=4, fig.width=6-----------------------------------------------
ggplot(data=NCline,aes(x=Time_from_present,y=Net_changes,color=Clusters))+
  geom_path(aes(x=Time_from_present,y=Net_changes,group=ID,color=Clusters),
            arrow = arrow(length = unit(0.10, "cm")))+
  facet_wrap(~sp_gender)+
  theme_classic()

## -----------------------------------------------------------------------------
Angl<-Angles
colnames(Angl)<-2:30
Angl$ID<-as.numeric(rownames(Angl))
Angl$cluster<-as.factor(groups)
Angl$sp_gender<-furseals$sp_gender[1:47]

Angline<- tidyr::pivot_longer(Angl, 1:29, 
                              names_to ="Time_from_present", 
                              values_to="Direction", 
                              names_transform = function(x) {as.numeric(x)-1}) 
colnames(Angline)[c(2,3)] = c("Clusters", "Group")

# range 15°
deg <- 15
# vector for range of direction of different bars
dir.breaks <- seq(0-(deg/2), 360+(deg/2), deg)
dir.binned <- cut(Angline$Direction,
                  breaks = dir.breaks,
                  ordered_result = TRUE)
# direction labels
dir.labels <- as.character(c(seq(0, 360-deg, by = deg), 0))
levels(dir.binned) <- dir.labels

# angles distribution in each range of direction
Angline$dir.binned <- dir.binned

# sort angles
df_sorted<-as.data.frame(table(Angline$dir.binned, Angline$Clusters))
colnames(df_sorted)<-c("dir.binned","Clusters","nb")
df_sorted = df_sorted[order(df_sorted$dir.binned),]

## ---- fig.height=4, fig.width=6-----------------------------------------------
ggplot(data=df_sorted, aes(x=dir.binned, y=nb, fill=Clusters)) +
  geom_bar(stat="identity")+
  scale_y_continuous(limits = c(0,110), expand = c(0, 0), 
                     breaks = c(0,25,50,75,110), 
                     labels = c(0,25,50,75,110)) +
  labs(x = 'Trajectory segment directions within fur seals clusters', y = 'number of trajectory segments') +
  coord_polar(start = -(deg/2)*(pi/180)) +
  theme_minimal()

## ----load Pike data-----------------------------------------------------------
data("pike")

## -----------------------------------------------------------------------------
Net_changes<-trajectoryLengths2D(pike[,7:8],pike$ID,pike$Time, relativeToInitial=TRUE) 
colnames(Net_changes)<-c("Net_changes", "Trajectory")
pike$Net_Changes<-Net_changes$Net_changes

## -----------------------------------------------------------------------------
D=dist(pike[,7:8])
Ds<-trajectoryDistances(D,pike$ID, surveys = NULL, distance.type = "DSPD",
                        symmetrization = "mean", add = TRUE, verbose = FALSE)

## ---- fig.height=4, fig.width=6-----------------------------------------------
Hst=3
colstd<-c("black","yellow","green","blue","grey","red")
hsxy <- hclust(Ds, "ward.D2")
plot(hsxy,hang = -1, main="distance Pike", cex=.6)
x<-rect.hclust (hsxy, h=Hst,
                border = colstd)
# Store clusters into a new data column
pike$Cluster<-cutree(hsxy, h=Hst)

## -----------------------------------------------------------------------------
Pike1<-pike[pike$Time %in% 1,]
Pike1<-Pike1[order(Pike1$ID, decreasing=FALSE),]
Pike1$Net_changes<-0
Pike2<-pike[pike$Time %in% 2,]
Pike2<-Pike2[order(Pike2$ID, decreasing=FALSE),]
Pike2$Net_changes<-Net_changes$Net_changes
data<-as.data.frame(rbind(Pike1,Pike2))

## ---- fig.height=4, fig.width=6-----------------------------------------------
ggplot(data=data,aes(x=d13C,y=d15N,shape=Trophic_status_initial))+
  geom_point(aes(size=Net_changes))+
  geom_path(aes(x=d13C,y=d15N,group=ID,color=factor(Cluster)),arrow = arrow(length = unit(0.30, "cm")))+
  geom_hline(yintercept=10, linetype="dashed", color = "black")+
  xlab(expression(delta^13*"C")) +
  ylab(expression(delta^15*"N"))+
  theme_minimal()

## ---- fig.height=3, fig.width=4-----------------------------------------------
gg_dist_d13C = ggplot(data, aes(d13C, fill=TimeL)) + geom_density(alpha=.5) 
gg_dist_d13C = gg_dist_d13C + ylab(expression(delta^13*"C"*" density"))
gg_dist_d13C = gg_dist_d13C + theme(axis.title.y=element_blank(),
                                    axis.text=element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks=element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.background =element_blank())
gg_dist_d13C = gg_dist_d13C + theme(legend.position = "none")
gg_dist_d13C + scale_x_continuous(limits = c(-33, -25))+scale_y_continuous(limits = c(0, 1))

## ---- fig.height=3, fig.width=4-----------------------------------------------
gg_dist_d15N = ggplot(data, aes(d15N, fill=TimeL)) + geom_density(alpha=.5) 
gg_dist_d15N = gg_dist_d15N + ylab(expression(delta^15*"N"*" density"))
gg_dist_d15N =gg_dist_d15N 
gg_dist_d15N =gg_dist_d15N + coord_flip()
gg_dist_d15N =gg_dist_d15N + theme(axis.title.y=element_blank(),
                                   axis.text=element_blank(),
                                   axis.line=element_blank(),
                                   axis.ticks=element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background =element_blank())
gg_dist_d15N =gg_dist_d15N +theme(legend.position = "none")
gg_dist_d15N + scale_x_continuous(limits = c(7, 14))+scale_y_continuous(limits = c(0, 1))

## ----load Alaska--------------------------------------------------------------
data("isoscape")

## -----------------------------------------------------------------------------
sites<-isoscape$station
surveys<-isoscape$Year
Angl<-trajectoryAngles2D(isoscape[,3:4],sites,surveys, betweenSegments = FALSE)
Length<-trajectoryLengths2D(isoscape[,3:4],sites,surveys)
data<-as.data.frame(cbind(isoscape[1:489,],Angl,Length[,1]))
colnames(data)<-c("Latitude","Longitude","d13C","d15N","Stations","Years","Angles","Lengths")
head(data)

## -----------------------------------------------------------------------------
angle<-data$Angles
Angles2<-c()
for (i in 1:length(angle)) {
  Angles2[i] <- c(ifelse(angle[i]==0,(angle[i]-270)*pi/180,
                         ifelse(angle[i]==180,(angle[i]-270)*pi/180,
                                ifelse(angle[i]==90,(angle[i]+270)*pi/180,
                                       ifelse(angle[i]==270,(angle[i]+270)*pi/180,
                                              ifelse(angle[i]==360,(angle[i]-270)*pi/180,  
                                                     ifelse(angle[i]>0 & angle[i]<90 ,(90-angle[i])*pi/180,
                                                            ifelse(angle[i]>90 & angle[i]<180 ,(90-angle[i])*pi/180,
                                                                   ifelse(angle[i]>180 & angle[i]<270,(180+(270-angle[i]))*pi/180,
                                                                          ifelse(angle[i]>270 & angle[i]<360,(90+(360-angle[i]))*pi/180,"ERROR"))))))))))
}

data$Angles2<-Angles2

## ---- fig.height=4, fig.width=6-----------------------------------------------
ggplot(data, 
          aes(x = Longitude, 
              y = Latitude, 
              fill = Lengths, 
              angle = Angles2, 
              radius = rescale(Lengths, c(0.3, 1)))) +
  geom_raster(interpolate = TRUE) +
  geom_spoke(arrow = arrow(length = unit(.07, 'inches'))) + 
  scale_fill_distiller(palette = "RdYlBu") + 
  coord_equal(expand = 0) + 
  theme(legend.position = 'bottom', 
        legend.direction = 'horizontal',
        panel.background = element_rect(fill = "white"))

## ----load heatmap, echo=T-----------------------------------------------------
data("heatmapdata")

## -----------------------------------------------------------------------------
head(heatmapdata)

## ---- echo=T------------------------------------------------------------------
#direction range
deg <- 15

dir.breaks <- c(0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360)


dir.binned <- cut(heatmapdata$Angles,
                  breaks = dir.breaks,
                  ordered_result = TRUE)

# bar labels
dir.labels <- as.character(c(seq(0, 360-deg, by = deg),0))

levels(dir.binned) <- dir.labels

heatmapdata$dir.binned <- dir.binned

data<-heatmapdata[,c(6,7,8,10)]

#direction vs SI patterns
data<-data[order(data$dir.binned, decreasing=FALSE),]
rownames(data)<-1:9206
data$ISpattern<- c(rep("+d13C/+d15N",2862),rep("+d13C/-d15N",1840),rep("-d13C/-d15N",2931), rep("-d13C/+d15N",1573))

data1<-as.data.frame(table(data$dir.binned,data$Years))

data2<-aggregate(x = data$Lengths, by = list(data$dir.binned, data$Years), FUN=sum, drop=FALSE)
data2[is.na(data2)] <- 0 


data1$Lengths<-data2$x
dfa<-data1
colnames(dfa)<-c("Directions","Periods","Nb_stations","Lengths")

## -----------------------------------------------------------------------------
head(dfa)

## ---- fig.height=4, fig.width=6-----------------------------------------------
ggplot(dfa, aes(Periods, Directions, fill= Nb_stations)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_minimal()+
  theme(axis.text.x = element_text(size=10, angle=90))

## ---- fig.height=3, fig.width=5-----------------------------------------------
df.Xbarplot<-aggregate(dfa$Lengths, by = list(dfa$Periods), FUN = sum)
colnames(df.Xbarplot)<-c("Periods","Lengths")
bp.x <- ggplot(data = df.Xbarplot, aes(x = factor(Periods), y = Lengths)) + 
  geom_bar(stat = "identity", aes(fill = Lengths)) + theme_minimal() +
  theme(axis.text.x = element_text(size = 10,angle=90), 
        axis.title.x = element_text(size = 20, margin = margin(10,0,0,0))) +
  labs(x = "Periods")
bp.x

## ---- fig.height=4, fig.width=4-----------------------------------------------
df.Ybarplot<-aggregate(dfa$Lengths, by = list(dfa$Directions), FUN = sum)
colnames(df.Ybarplot)<-c("Directions","Lengths")
df.Ybarplot$ISpattern<- c(rep("+d13C/+d15N",6),rep("+d13C/-d15N",6),rep("-d13C/-d15N",6), rep("-d13C/+d15N",6))


bp.y <- ggplot(data = df.Ybarplot, aes(x = factor(Directions), y = Lengths,fill = ISpattern)) + 
  geom_bar(stat="identity") + theme_minimal() + coord_flip()
bp.y

