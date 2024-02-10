###Open necessary packages####
library(vegan)
library(tidyr)
library(dplyr)
library(readr)
library(plyr)
library(tidyverse)
library(Rmisc)
library(ggplot2)
library(emmeans)
library(nlme)
library(lattice)
library(multcompView)
library(multcomp)
library(clipr)
library(MuMIn)

###Read in csv of IC Agronomic Data and Organize####
LEG<-read.csv("D:/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/LegDat_R_V11_JJ.csv")
LEG<-read.csv("/Volumes/GoogleDrive/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/LegDat_R_V11_JJ.csv")
#LEG$Trt_code<-revalue(LEG$Trt_code, c("AlCl"="Alsike Clover", "BTre"="Birdsfoot Trefoil", "Alf"="Alfalfa", "CMV"="Canada Milk Vetch", "HiN"="90 kg ha", "MedN"="45 kg ha","None"="Control", "RdCl"="Red Clover", "WhCl"="White Clover"))
#LEG$Trt_code<-factor(LEG$Trt_code, levels=c("Control", "45 kg ha", "90 kg ha", "Alsike Clover", "Birdsfoot Trefoil", "Alfalfa", "Canada Milk Vetch", "Red Clover", "White Clover"))
LEG$Trt_code<-revalue(LEG$Trt_code, c("AlCl"="A. clover", "BTre"="B. trefoil", "Alf"="Alfalfa", "CMV"="C. milkvetch", "HiN"="90", "MedN"="45","None"="Control", "RdCl"="R. clover", "WhCl"="W. clover"))
LEG$Trt_code<-factor(LEG$Trt_code, levels=c("Control", "45", "90", "A. clover", "B. trefoil", "Alfalfa", "C. milkvetch", "R. clover", "W. clover"))


#Ensure factors as treated as such
LEG$fYear<-as.factor(LEG$Year)
LEG$fRep<-as.factor(LEG$Rep)
LEG$fPlot<-as.factor(LEG$Plot)
LEG$Grain.kg.ha <- LEG$Recalculated.IWG.Grain.kg.ha
LEG$Straw.kg.ha <- LEG$Recalculated.IWG.Straw.kg.ha
LEG$Leg.bio.kg.ha <- LEG$Legume.Biomass.kg.ha
LEG$Weed.bio.kg.ha <-LEG$Weed.Biomass.kg.ha
LEG$LegWeed.bio.kg.ha <-LEG$Legume.weed.Biomass.kg.ha
LEG$NTrans <- LEG$X..N_transfer
LEG$fStandAge <- as.factor(LEG$Year-2016)
LEG$StandAge <- LEG$Year-2016
LEG$LRR<-LEG$Leg.bio.kg.ha/(LEG$Straw.kg.ha)#+LEG$Grain.kg.ha)

SPdat<-subset(LEG, Location=="Saint Paul")
Lamdat<-subset(LEG, Location=="Lamberton")
Rosdat<-subset(LEG, Location=="Rosemount")

###Summary SE to find actual means####
realmean<-summarySE(subset(LEG, Trt_code!="X"), "Grain.kg.ha", c("fYear", "Trt_code", "Location"), na.rm=TRUE)

#Legume biomass / grain yield ####

sumLRR2<-summarySE(subset(LEG, fYear=="2017"), "Leg.bio.kg.ha", c("Location", "fRep", "Trt_code"), na.rm=TRUE)
sumLRR2<-rbind(sumLRR2, summarySE(subset(LEG, fYear=="2017"), "Leg.bio.kg.ha", c("Location", "fRep","Trt_code"), na.rm=TRUE))
sumLRR2<-rbind(sumLRR2, summarySE(subset(LEG, fYear=="2017"), "Leg.bio.kg.ha", c("Location","fRep", "Trt_code"), na.rm=TRUE))
sumLRR2$fYear<-as.factor(rep(c("2017","2018","2019"), each=106))
sumLRR2$Grain.kg.ha<-summarySE(LEG, "Grain.kg.ha", c("fYear", "Location","fRep", "Trt_code"), na.rm=TRUE)[,6]
sumLRR2$NTrans<-summarySE(LEG, "NTrans", c("fYear", "Location","fRep", "Trt_code"), na.rm=TRUE)[,6]
sumLRR2$IWG_CN<-summarySE(LEG, "IWG_CN", c("fYear", "Location","fRep", "Trt_code"), na.rm=TRUE)[,6]

#Correlation between 2017 legume biomass and grain yield in 2017, 2018, and 2019
regmod1<-lme(Grain.kg.ha~Leg.bio.kg.ha*fYear,random = ~1|Location/fRep, na.action=na.omit, data=sumLRR2)
anova(regmod1)
summary(regmod1)
plot(regmod1)
r.squaredGLMM(regmod1)
summary(lme(Grain.kg.ha~Leg.bio.kg.ha,random = ~1|Location/fRep, na.action=na.omit,
            data=subset(sumLRR2, fYear=="2018")))

ggplot(data=sumLRR2, aes(x=Leg.bio.kg.ha, y=Grain.kg.ha)) + 
  geom_point(size=1.5)+
  facet_grid(~fYear) +
  geom_smooth(se=F, method="lm", color="black")+
  ylim(c(0,900))+
  ylab(expression("Grain yield " ~ (kg ~ ha^{-1})))+ 
  xlab(expression("Legume biomass " ~ (kg ~ ha^{-1})))+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.text.x=element_text(size=10, color='black'))
ggsave("FigureX.png", width=16.5, height=8.5, units="cm", path="D:/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")
ggsave("FigureX.png", width=16.5, height=8.5, units="cm", path="/Volumes/GoogleDrive/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

#Correlation between 2017 legume biomass and n transfer in 2017, 2018, and 2019

regmod2<-lme(NTrans~Leg.bio.kg.ha*fYear,random = ~1|Location/fRep, na.action=na.omit, data=sumLRR2)
anova(regmod2)
summary(regmod2)
plot(regmod2)
r.squaredGLMM(regmod2)

regmod2<-lme(NTrans~Leg.bio.kg.ha*fYear,random = ~1|Location/fRep, na.action=na.omit, data=sumLRR2)
anova(regmod2)
summary(regmod2)
ggplot(data=sumLRR2, aes(x=Leg.bio.kg.ha, y=NTrans)) + 
  geom_point(size=1.5)+
  facet_grid(~fYear) +
  geom_smooth(se=F, method="lm", color="black")+
  #ylim(c(0,900))+
  ylab("N transfer")+ 
  xlab(expression("Legume biomass " ~ (kg ~ ha^{-1})))+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.text.x=element_text(size=10, color='black'))

#Look at the slope for each year
summary(lme(NTrans~Leg.bio.kg.ha,random = ~1|Location/fRep, 
            na.action=na.omit, data=subset(sumLRR2, fYear=="2019")))

#Correlation between 2017 legume biomass and IWG C:N in 2017, 2018, and 2019
regmod3<-lme(IWG_CN~Leg.bio.kg.ha*fYear,random = ~1|Location/fRep, na.action=na.omit, data=sumLRR2)
anova(regmod3)
summary(regmod3)


#Correlation between 2017 legume biomass and 2019 N transfer
summary(lme(NTrans~Leg.bio.kg.ha,random = ~1|Location/fRep, na.action=na.omit,
            data=subset(sumLRR2, fYear=="2019")))

#Correlation between N transfer and grain yield
regmod4<-lme(Grain.kg.ha~NTrans*fYear,random = ~1|Location/fRep, na.action=na.omit, data=LEG)
anova(regmod4)
summary(regmod4)
r.squaredGLMM(regmod4)

ggplot(data=LEG, aes(y=Grain.kg.ha, x=NTrans)) + 
  geom_point(size=1.5)+
  facet_grid(~fYear) +
  geom_smooth(se=F, method="lm", color="black")+
  ylim(c(0,900))+
  ylab(expression("Grain yield " ~ (kg ~ ha^{-1})))+ 
  xlab("N Transfer")+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.text.x=element_text(size=10, color='black'))

####Within-year correlations ####
#2017
cordat<-subset(LEG, Year=="2017")[,c(28,29,30,12,31,21,32)]
names(cordat)<-c("IWG grain","IWG biomass","Legume biomass", "Total biomass", "Weed biomass", "IWG C:N", "N transfer")
pcor<-round(cor(cordat, use = "na.or.complete"),2)
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(pcor)
library(reshape2)
mpcor<-melt(upper_tri)
library(ggplot2)
ggplot(data = mpcor, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "green", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))+
  coord_fixed()+
  ggtitle("A")+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  theme(
    plot.title=element_text(size=14,face='bold', hjust=.02),
    axis.title.x = element_blank(),
    axis.text=element_text(color="black", size=14),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    legend.background = element_rect(color="white"),
    legend.text=element_text(size=14),
    legend.position=c(0.2,0.8),
    legend.title=element_text(size=14),
    axis.ticks = element_blank())
ggsave("Leg2017.png", width=7, height=6, units="in", path="/Volumes/GoogleDrive/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

#2018
cordat<-subset(LEG, Year=="2018")[,c(28,29,30,12,31,21,32)]
names(cordat)<-c("IWG grain","IWG biomass","Legume biomass", "Total biomass", "Weed biomass", "IWG C:N", "N transfer")
pcor<-round(cor(cordat, use = "na.or.complete"),2)
upper_tri <- get_upper_tri(pcor)
mpcor<-melt(upper_tri)
ggplot(data = mpcor, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "green", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))+
  coord_fixed()+
  ggtitle("B")+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  theme(
    plot.title=element_text(size=14,face='bold', hjust=.02),
    axis.title.x = element_blank(),
    axis.text.y=element_blank(),
    axis.text=element_text(color="black", size=14),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    legend.position="none",
    axis.ticks = element_blank())
ggsave("Leg2018.png", width=7, height=6, units="in", path="/Volumes/GoogleDrive/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

#2019
#201
cordat<-subset(LEG, Year=="2019")[,c(28,29,30,12,31,21,32)]
names(cordat)<-c("IWG grain","IWG biomass","Legume biomass", "Total biomass", "Weed biomass", "IWG C:N", "N transfer")
pcor<-round(cor(cordat, use = "na.or.complete"),2)
upper_tri <- get_upper_tri(pcor)
mpcor<-melt(upper_tri)
ggplot(data = mpcor, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "green", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))+
  coord_fixed()+
  ggtitle("C")+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  theme(
    plot.title=element_text(size=14,face='bold', hjust=.02),
    axis.title.x = element_blank(),
    axis.text.y=element_blank(),
    axis.text=element_text(color="black", size=14),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    legend.position="none",
    axis.ticks = element_blank())
ggsave("Leg2019.png", width=7, height=6, units="in", path="/Volumes/GoogleDrive/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

###Total Biomass Analysis####
histogram(~Total.VBiomass|Location, data=LEG)
TotVBioMod1<-lme(Total.VBiomass~Trt_code*fStandAge,random = ~1|fRep, na.action=na.omit, data=SPdat)
anova(TotVBioMod1)
cld(emmeans(TotVBioMod1, ~fStandAge))
cld(emmeans(TotVBioMod1, ~Trt_code|fStandAge))
cld(emmeans(TotVBioMod1, ~Trt_code))

TotVBioMod2<-lme(Total.VBiomass~Trt_code*fStandAge,random = ~1|fRep, na.action=na.omit, data=Lamdat)
anova(TotVBioMod2)
summary(TotVBioMod2)
cld(emmeans(TotVBioMod2, ~fStandAge))
cld(emmeans(TotVBioMod2, ~Trt_code|fStandAge))
cld(emmeans(TotVBioMod2, ~Trt_code))

TotVBioMod3<-lme(Total.VBiomass~Trt_code*fStandAge,random = ~1|fRep, na.action=na.omit, data=Rosdat)
anova(TotVBioMod3)
cld(emmeans(TotVBioMod3, ~fStandAge))
cld(emmeans(TotVBioMod3, ~Trt_code|fStandAge))
cld(emmeans(TotVBioMod3, ~Trt_code), order=F)

### Weed Biomass Response ####
#SP
histogram(~sqrt(Weed.bio.kg.ha)|Location, data=LEG)
WeedBioMod1<-lme(sqrt(Weed.bio.kg.ha)~Trt_code*fStandAge, random = ~1|fRep, na.action=na.omit, data=SPdat)
anova(WeedBioMod1)
cld(emmeans(lm(Weed.bio.kg.ha~Trt_code*fStandAge, data=SPdat), ~Trt_code))
cld(emmeans(lm(Weed.bio.kg.ha~Trt_code*fStandAge, data=SPdat), ~fStandAge))


#Lamberton #Use sqrt response for mean comparison letters, but no transformation for mean estimate
WeedBioMod2<-lme(sqrt(Weed.bio.kg.ha)~Trt_code*fStandAge, random = ~1|fRep, na.action=na.omit, data=Lamdat)
anova(WeedBioMod2)
cld(emmeans(lme(Weed.bio.kg.ha~Trt_code*fStandAge, random = ~1|fRep, na.action=na.omit, data=Lamdat), ~fStandAge))
cld(emmeans(lme(sqrt(Weed.bio.kg.ha)~Trt_code*fStandAge, random = ~1|fRep, na.action=na.omit, data=Lamdat), ~fStandAge))

cld(emmeans(lm(Weed.bio.kg.ha~Trt_code*fStandAge, data=Lamdat), ~fStandAge))

#Rosemount
WeedBioMod3<-lme(sqrt(Weed.bio.kg.ha)~Trt_code*fStandAge, random = ~1|fRep, na.action=na.omit, data=Rosdat)
WeedBioMod4<-lm(Weed.bio.kg.ha~Trt_code*fStandAge, data=Rosdat)
anova(WeedBioMod3)
cld(emmeans(WeedBioMod4, ~Trt_code))
cld(emmeans(WeedBioMod4, ~fStandAge))


###Grain Yield####
#Rosemount
#Stand Age as Factor
RosGrainmod2<-lme(Grain.kg.ha~Trt_code*fStandAge, data=Rosdat, random = ~1|fRep, na.action = na.omit)
anova(RosGrainmod2)
cld(emmeans(RosGrainmod2, ~fStandAge))
cld(emmeans(RosGrainmod2, ~Trt_code|fStandAge))
cld(emmeans(RosGrainmod2, ~Trt_code))

glht(RosGrainmod2, linfct = mcp(Trt_code = c("45 - Control= 0")))
glht(RosGrainmod2, linfct = mcp(Trt_code = "Tukey"))

RosGrainmod3<-lm(Grain.kg.ha~Trt_code*fStandAge*Leg.bio.kg.ha, data=Rosdat)
d<-anova(RosGrainmod3)
summary(RosGrainmod3)

#Lamberton
#Stand Age as Factor
LamGrainmod2<-lme(Grain.kg.ha~Trt_code*fStandAge, data=Lamdat, random = ~1|fRep, na.action = na.omit)
anova(LamGrainmod2)
cld(emmeans(LamGrainmod2, ~fStandAge))
cld(emmeans(LamGrainmod2, ~Trt_code|fStandAge))
cld(emmeans(LamGrainmod2, ~Trt_code))

#Saint Paul
#Stand Age as Factor
SPGrainmod2<-lme(Grain.kg.ha~Trt_code*fStandAge, data=SPdat, random = ~1|fRep, na.action = na.omit)
anova(SPGrainmod2)
cld(emmeans(SPGrainmod2, ~fStandAge))
cld(emmeans(SPGrainmod2, ~Trt_code|fStandAge))
cld(emmeans(SPGrainmod2, ~Trt_code))

#Lamberton grain yields Figure 1 ####
sumLEGL<-summarySE(Lamdat, "Grain.kg.ha", c("fYear", "Trt_code", "fStandAge"), na.rm=TRUE)
sumLEGL2<-merge(sumLEGL,
                      as.data.frame(cld(emmeans(LamGrainmod2, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                      by=c("Trt_code","fStandAge"))
sumLEGL2$.group[c(which(sumLEGL2$fStandAge=="2"))]<-NA
pps = position_dodge(width = .75)

ggplot(data=sumLEGL2, aes(x=Trt_code, y=Grain.kg.ha)) + 
  geom_bar(size=1.5, stat="identity", position = pps) + 
  facet_grid(~fYear) +
  geom_errorbar(aes(ymin=Grain.kg.ha-se, ymax=Grain.kg.ha+se), width=0.2, position = pps) +
  xlab("Treatment") + 
  ylab(expression("Grain yield " ~ (kg ~ ha^{-1})))+ 
  ylim(0,750) +
  geom_text(aes(x=Trt_code, y=Grain.kg.ha+se+25, label=.group), 
            show.legend = F, size=2.75, color="black", position = position_dodge(1),
            hjust="center")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=10, angle=45,
                                 hjust=0.95, color='black'))
ggsave("Figure1.png", width=16.5, height=8.5, units="cm", path="D:/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

###IWG Biomass####
#Stand Age as factor
RosStrawmod3<-lme(Straw.kg.ha~Trt_code*fStandAge, data=Rosdat, random = ~1|fRep, na.action = na.omit)
anova(RosStrawmod3)
cld(emmeans(RosStrawmod3, ~fStandAge))
cld(emmeans(RosStrawmod3, ~Trt_code))
cld(emmeans(RosStrawmod3, ~Trt_code|fStandAge))


#Lamberton
#Stand Age as factor
LamStrawmod3<-lme(Straw.kg.ha~Trt_code*fStandAge, data=Lamdat, random = ~1|fRep, na.action = na.omit)
anova(LamStrawmod3)
cld(emmeans(LamStrawmod3, ~fStandAge))
cld(emmeans(LamStrawmod3, ~Trt_code))
cld(emmeans(LamStrawmod3, ~Trt_code|fStandAge))

#Saint Paul
#Stand Age as factor
SPStrawmod4<-lme(Straw.kg.ha~Trt_code*fStandAge, data=SPdat, random = ~1|fRep, na.action = na.omit)
anova(SPStrawmod4)
cld(emmeans(SPStrawmod4, ~fStandAge))
cld(emmeans(SPStrawmod4, ~Trt_code))
cld(emmeans(SPStrawmod4, ~Trt_code|fStandAge))



###Legume Biomass Yield####
#Rosemount
RosLBiomod1<-lme(Leg.bio.kg.ha~Trt_code*fStandAge, data=Rosdat, random = ~1|fRep, na.action = na.omit)
anova(RosLBiomod1)
cld(emmeans(RosLBiomod1, ~Trt_code))
cld(emmeans(RosLBiomod1, ~fStandAge))
#cld(lstrends(RosLBiomod1,~Trt_code, var="StandAge"))

#Lamberton
LamLBiomod1<-lme(Leg.bio.kg.ha~Trt_code*fStandAge, data=Lamdat, random = ~1|fRep, na.action = na.omit)
anova(LamLBiomod1)
cld(emmeans(LamLBiomod1, ~Trt_code))
cld(emmeans(LamLBiomod1, ~fStandAge))
#Saint Paul
SPLBiomod1<-lme(Leg.bio.kg.ha~Trt_code*fStandAge, data=SPdat, random = ~1|fRep, na.action = na.omit)
anova(SPLBiomod1)
cld(emmeans(SPLBiomod1, ~Trt_code))
cld(lstrends(SPLBiomod1,~Trt_code, var="StandAge"))

#Legume Yield Summary
sumLBIO<-summarySE(subset(LEG, Trt_code!="X"&Legume=="1"), "Leg.bio.kg.ha", c("fYear", "Trt_code", "Location", "fStandAge"), na.rm=TRUE)
sumLBIO$Leg.bio.kg.ha[which(sumLBIO$Leg.bio.kg.ha=="NaN")]<-NA
sumLBIO_letters<-rbind(as.data.frame(cld(emmeans(RosLBiomod1, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                       as.data.frame(cld(emmeans(LamLBiomod1, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                       as.data.frame(cld(emmeans(SPLBiomod1, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)))
sumLBIO_letters$Location<-rep(c("Rosemount","Lamberton","Saint Paul"), each=18)
sumLBIO_letters$.group[37:54]<-NA
sumLBIO_letters$.group[c(which(sumLBIO_letters$fStandAge=="1"&sumLBIO_letters$Location=="Rosemount"))]<-NA

sumLBIO2<-merge(sumLBIO,
                sumLBIO_letters,
                by=c("Trt_code","fStandAge","Location"))
#Figure
ggplot(data=sumLBIO2, aes(x=Trt_code, y=Leg.bio.kg.ha)) + 
  geom_bar(size=1.5, stat="identity", position = pps) + 
  facet_grid(Location~fYear) +
  geom_errorbar(aes(ymin=Leg.bio.kg.ha-se, ymax=Leg.bio.kg.ha+se), width=0.2, position = pps) +
  xlab("Treatment") + 
  ylab(expression("Legume biomass yield " ~ (kg ~ ha^{-1})))+ 
  #ylim(0,750) +
  geom_text(aes(x=Trt_code, y=Leg.bio.kg.ha+se+300, label=.group), 
            show.legend = F, size=2.75, color="black", position = position_dodge(1),
            hjust="center")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=10, angle=45,
                                 hjust=0.95, color='black'))
ggsave("Figure2.png", width=16.5, height=15, units="cm", path="D:/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

#Isotope Analyses


#Stand Age As Factor
f15N1<-lm(IWG_15N~Trt_code*fStandAge, data=SPdat, na.action = na.omit)
anova(f15N1)
cld(emmeans(f15N1, ~Trt_code))
cld(emmeans(f15N1, ~fStandAge))
cld(emmeans(f15N1, ~Trt_code|fStandAge))

f15N2<-lme(IWG_15N~Trt_code*fStandAge, data=Lamdat, random = ~1|fRep, na.action = na.omit)
anova(f15N2)
cld(emmeans(f15N2, ~Trt_code))
cld(emmeans(f15N2, ~fStandAge))
cld(emmeans(f15N2, ~Trt_code|fStandAge))

#Stand Age As Factor SP
f15N3<-lme(IWG_15N~Trt_code*fStandAge, data=Rosdat, random = ~1|fRep,na.action = na.omit)
anova(f15N3)
cld(emmeans(f15N3, ~Trt_code))
cld(emmeans(f15N3, ~fStandAge))


#15N Ratio Summary and Figure
sumN<-summarySE(subset(LEG, Trt_code!="X"), "IWG_15N", c("fYear", "Trt_code", "Location", "fStandAge"), na.rm=TRUE)
sumN$IWG_15N[which(sumN$IWG_15N=="NaN")]<-NA

sumN_letters<-rbind(as.data.frame(cld(emmeans(f15N3, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                       as.data.frame(cld(emmeans(f15N2, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                       as.data.frame(cld(emmeans(f15N1, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)))
sumN_letters$Location<-rep(c("Rosemount","Lamberton","Saint Paul"), each=27)
sumN_letters$.group[1:27]<-NA
sumN_letters$.group[c(which(sumN_letters$fStandAge=="3"&sumN_letters$Location=="Lamberton"))]<-NA
sumN_letters$.group[c(which(sumN_letters$fStandAge=="3"&sumN_letters$Location=="Saint Paul"))]<-NA
sumN_letters$.group[c(which(sumN_letters$fStandAge=="2"&sumN_letters$Location=="Saint Paul"))]<-NA

sumN2<-merge(sumN,
             sumN_letters,
                by=c("Trt_code","fStandAge","Location"))

ggplot(data=sumN2, aes(x=Trt_code, y=IWG_15N)) + 
  geom_bar(size=1.5, stat="identity", position = pps) + 
  facet_grid(Location~fYear) +
  geom_errorbar(aes(ymin=IWG_15N-se, ymax=IWG_15N+se), width=0.2, position = pps) +
  xlab("Treatment") + 
  ylab(expression(paste("IWG ",delta^{15}, "N (\u2030)")))+ 
  #ylim(0,750) +
  geom_text(aes(x=Trt_code, y=IWG_15N+se+.7, label=.group), 
            show.legend = F, size=2.75, color="black", position = position_dodge(1),
            hjust="center")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=10, angle=45,
                                 hjust=0.95, color='black'))
ggsave("Figure4.png", width=16.5, height=15, units="cm", path="D:/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

###Carbon/Nitrogen Ratio####
#Rosemount
Nmod3<-lme(IWG_CN~Trt_code*fStandAge, data = Rosdat, random = ~1|fRep, na.action = na.omit)
anova(Nmod3)
cld(emmeans(Nmod3, ~Trt_code))
cld(lstrends(Nmod3,~fStandAge))

#Lamberton
Nmod2<-lme(IWG_CN ~ Trt_code * fStandAge, data = Lamdat, random = ~1|fRep, na.action = na.omit)
anova(Nmod2)
cld(emmeans(Nmod2, ~Trt_code))
cld(lstrends(Nmod2,~fStandAge)) 
  
#Saint Paul
Nmod4<-lm(IWG_CN ~ Trt_code * fStandAge, data = SPdat, na.action = na.omit)
anova(Nmod4)
cld(emmeans(Nmod4, ~Trt_code))
cld(emmeans(Nmod4,~fStandAge)) 


#Lamberton C:N Figure 3
sumLCN<-summarySE(Lamdat, "IWG_CN", c("fYear", "Trt_code", "fStandAge"), na.rm=TRUE)
sumLCN2<-merge(sumLCN,
                as.data.frame(cld(emmeans(Nmod2, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                by=c("Trt_code","fStandAge"))
sumLCN2$.group[c(which(sumLCN2$fStandAge=="1"))]<-NA


ggplot(data=sumLCN2, aes(x=Trt_code, y=IWG_CN)) + 
  geom_bar(size=1.5, stat="identity", position = pps) + 
  facet_grid(~fYear) +
  geom_errorbar(aes(ymin=IWG_CN-se, ymax=IWG_CN+se), width=0.2, position = pps) +
  xlab("Treatment") + 
  ylab("IWG C:N")+ 
  geom_text(aes(x=Trt_code, y=IWG_CN+se+5, label=.group), 
            show.legend = F, size=2.75, color="black", position = position_dodge(1),
            hjust="center")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=5),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=10, angle=45,
                                 hjust=0.95, color='black'))
ggsave("Figure3.png", width=16.5, height=8.5, units="cm", path="D:/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

###N Transfer####
#Rosemount
NTmod3<-lme(NTrans ~ Trt_code * fStandAge, data = subset(LEG, Location=="Rosemount"&Legume=="1"), random = ~1|fRep, na.action = na.omit)
anova(NTmod3)
cld(emmeans(NTmod3, ~Trt_code)) 

NTmodR1<-lme(NTrans ~ (Trt_code-1) * fStandAge, data = subset(LEG, Location=="Rosemount"&Legume=="1"), random = ~1|fRep, na.action = na.omit)
anova(NTmodR1)
summary(NTmodR1)

#Lamberton
NTmod2<-lme(NTrans ~ Trt_code * fStandAge, data = subset(LEG, Location=="Lamberton"&Legume=="1"), random = ~1|fRep, na.action = na.omit)
anova(NTmod2)
cld(emmeans(NTmod2, ~Trt_code))
cld(emmeans(NTmod2, ~fStandAge))

NTmodL1<-lme(NTrans ~ (Trt_code-1) * fStandAge, data = subset(LEG, Location=="Lamberton"&Legume=="1"), random = ~1|fRep, na.action = na.omit)
anova(NTmodL1)
summary(NTmodL1)

anova(lme(NTrans ~ (Trt_code-1), data = subset(LEG, fYear=="2017"&Location=="Lamberton"&Legume=="1"), random = ~1|fRep, na.action = na.omit))
summary(lme(NTrans ~ (Trt_code-1), data = subset(LEG, fYear=="2017"&Location=="Lamberton"&Legume=="1"), random = ~1|fRep, na.action = na.omit))

#Saint Paul
NTmod4<-lm(NTrans ~ Trt_code * fStandAge, data = subset(LEG, Location=="Saint Paul"&Legume=="1"), na.action = na.omit)
anova(NTmod4)
cld(emmeans(NTmod4, ~Trt_code))
cld(emmeans(NTmod4, ~fStandAge))

NTmodSP1<-lm(NTrans ~ (Trt_code-1) * fStandAge, data = subset(LEG, Location=="Saint Paul"&Legume=="1"), na.action = na.omit)
anova(NTmodSP1)
summary(NTmodSP1)

summary(lme(NTrans ~ (Trt_code-1), data = subset(LEG, fYear=="2019"&Location=="Saint Paul"&Legume=="1"), random = ~1|fRep, na.action = na.omit))

#N Transfer Summary and Figure
sumNT<-summarySE(subset(LEG, Trt_code!="X"&Legume=="1"), "NTrans", c("fYear", "Trt_code", "Location", "fStandAge"), na.rm=TRUE)
sumNT$NTrans[which(sumNT$NTrans=="NaN")]<-NA

sumNT_letters<-rbind(as.data.frame(cld(emmeans(NTmod3, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                    as.data.frame(cld(emmeans(NTmod2, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)),
                    as.data.frame(cld(emmeans(NTmod4, ~Trt_code|fStandAge), Letters=LETTERS, reversed=T)))
sumNT_letters$Location<-rep(c("Rosemount","Lamberton","Saint Paul"), each=18)
sumNT_letters$.group[1:18]<-NA
sumNT_letters$.group[c(which(sumNT_letters$fStandAge=="2"&sumNT_letters$Location=="Saint Paul"))]<-NA
sumNT_letters$.group[c(which(sumNT_letters$fStandAge=="3"&sumNT_letters$Location=="Saint Paul"))]<-NA
sumNT_letters$.group[c(which(sumNT_letters$fStandAge=="3"&sumNT_letters$Location=="Lamberton"))]<-NA

sumNT2<-merge(sumNT,
             sumNT_letters,
             by=c("Trt_code","fStandAge","Location"))

ggplot(data=sumNT2, aes(x=Trt_code, y=NTrans*100)) + 
  geom_point(size=1.5)+#, stat="identity", position = pps) + 
  facet_grid(Location~fYear) +
  geom_hline(yintercept = 0, color="black")+
  geom_errorbar(aes(ymin=NTrans*100-se*100, ymax=NTrans*100+se*100), width=0.2, position = pps) +
  xlab("Treatment") +
  ylab("N transfer (%)")+
  #ylim(0,750) +
  geom_text(aes(x=Trt_code, y=NTrans*100+se*100+20, label=.group), 
            show.legend = F, size=2.75, color="black", position = position_dodge(1),
            hjust="center")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=10, angle=45,
                                 hjust=0.95, color='black'))
ggsave("Figure5.png", width=16.5, height=15, units="cm", path="D:/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")

#Subset of full figure
sumNT3<-sumNT2
#manually indicate which means are different from 0. Use models above
sumNT3$zdif<-NA
sumNT3$zdif[c(4,7,12,25,39,40,43)]<-"*"
sumNT3$Env<-paste0(sumNT3$Location,":",sumNT3$fYear)
sumNT3<-subset(sumNT3, Env=="Saint Paul:2017"|Env=="Lamberton:2018"|Env=="Lamberton:2019")
sumNT3$Env<-factor(sumNT3$Env, levels=c("Saint Paul:2017", "Lamberton:2018", "Lamberton:2019"))

ggplot(data=sumNT3, aes(x=Trt_code, y=NTrans*100)) + 
  geom_point(size=1.5)+#, stat="identity", position = pps) + 
  facet_grid(~Env) +
  geom_hline(yintercept = 0, color="black")+
  geom_errorbar(aes(ymin=NTrans*100-se*100, ymax=NTrans*100+se*100), width=0.2, position = pps) +
  xlab("Treatment") +
  ylab("N transfer (%)")+
  #ylim(0,750) +
  coord_cartesian(ylim=c(-75,75))+
  scale_y_continuous(breaks=c(-75, -50, -25, 0, 25, 50, 75))+
  geom_text(aes(x=Trt_code, y=NTrans*100+se*100+10, label=zdif), 
            show.legend = F, size=7, color="black", position = position_dodge(1),
            hjust="center")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=10),
        strip.background =element_rect(color="black", fill="white"),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=10, angle=45,
                                 hjust=0.95, color='black'))
ggsave("Figure5v2.png", width=16.5, height=8.5, units="cm", path="/Volumes/GoogleDrive/My Drive/Mac/Students/Evelyn_Reiley/Legume_MS/")


###Difference in 15N IWG - 15N Legume####
#Rosemount
Difmod3<-lme(Diff_IWG_LEG ~ Trt_code * StandAge, data = subset(LEG, Location=="Rosemount"), random = ~1|fRep, na.action = na.omit)
anova(Difmod3)
summary(Difmod3)
cld(emmeans(Difmod3, ~Trt_code))

#Lamberton
Difmod2<-lme(Diff_IWG_LEG ~ Trt_code * StandAge, data = subset(LEG, Location=="Lamberton"), random = ~1|fRep, na.action = na.omit)
anova(Difmod2)
summary(Difmod2)
cld(emmeans(Difmod2, ~Trt_code))

#Saint Paul
Difmod4<-lme(Diff_IWG_LEG ~ Trt_code * StandAge, data = subset(LEG, Location=="Saint Paul"), random = ~1|fRep, na.action = na.omit)
anova(Difmod4)
summary(Difmod4)
cld(emmeans(Difmod4, ~Trt_code))

#All Sites
Difmod1<-lme(Diff_IWG_LEG ~ Trt_code * StandAge, data = LEG, random = ~1|fRep/Location, na.action = na.omit)
anova(Difmod1)
summary(Difmod1)
cld(emmeans(Difmod1, ~Trt_code))

#15N Difference Summary and Figure
sumDiff<-summarySE(subset(LEG, Trt_code!="X"), "Diff_IWG_LEG", c("fYear", "Trt_code", "Location"), na.rm=TRUE)
sumDiff$Diff_IWG_LEG[which(sumDiff$Diff_IWG_LEG=="NaN")]<-NA

ggplot(data=sumDiff, aes(x=Trt_code, y=Diff_IWG_LEG, fill=Trt_code)) + 
  geom_bar(size=1.5, stat="identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin=Diff_IWG_LEG-se, ymax=Diff_IWG_LEG+se)) +
  facet_grid(Location~fYear) +
  ylab("15N IWG - LEG Difference") + 
  #ylim(-5,12) +
  ggtitle("Legume Trial IWG-LEG 15N Difference") +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(plot.title=element_text(size=15,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=12),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.title=element_blank(),
        legend.key = element_rect(colour = "transparent", fill = NA),
        legend.text=element_text(size=12),
        legend.key.size = unit(6,"mm"))

###Fraction of N from Legume####
#Rosemount
LFmod3<-lme(Leg_Fraction ~ Trt_code * StandAge, data = subset(LEG, Location=="Rosemount"), random = ~1|fRep, na.action = na.omit)
anova(LFmod3)
summary(LFmod3)
cld(emmeans(LFmod3, ~Trt_code))

#Lamberton
LFmod2<-lme(Leg_Fraction ~ Trt_code * StandAge, data = subset(LEG, Location=="Lamberton"), random = ~1|fRep, na.action = na.omit)
anova(LFmod2)
summary(LFmod2)
cld(emmeans(LFmod2, ~Trt_code))

#Saint Paul
LFmod4<-lme(Leg_Fraction ~ Trt_code * StandAge, data = subset(LEG, Location=="Saint Paul"), random = ~1|fRep, na.action = na.omit)
anova(LFmod4)
summary(LFmod4)
cld(emmeans(LFmod4, ~Trt_code))

#All Sites
LFmod1<-lme(Leg_Fraction ~ Trt_code * StandAge, data = LEG, random = ~1|fRep/Location, na.action = na.omit)
anova(LFmod1)
summary(LFmod1)
cld(emmeans(LFmod1, ~Trt_code))

#Legume N Fraction Summary and Figure
sumLF<-summarySE(subset(LEG, Trt_code!="X"), "Leg_Fraction", c("fYear", "Trt_code", "Location"), na.rm=TRUE)
sumLF$Leg_Fraction[which(sumLF$Leg_Fraction=="NaN")]<-NA

ggplot(data=sumLF, aes(x=Trt_code, y=Leg_Fraction, fill=Trt_code)) + 
  geom_bar(size=1.5, stat="identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin=Leg_Fraction-se, ymax=Leg_Fraction+se)) +
  facet_grid(Location~fYear) +
  ylab("Fraction Legume Derived") + 
  #ylim(-5,12) +
  ggtitle("IWG Legume-derived N Fraction") +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(plot.title=element_text(size=15,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text.y = element_text(size=12),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        legend.title=element_blank(),
        legend.key = element_rect(colour = "transparent", fill = NA),
        legend.text=element_text(size=12),
        legend.key.size = unit(6,"mm"))



###Correlations####
sumCOR<-summarySE(subset(LEG, Trt_code!="X"), "Grain.kg.ha", c("fYear", "Trt_code", "Location", "Leg_Fraction", "Diff_IWG_LEG", "Straw.kg.ha", "Leg.bio.kg.ha"), na.rm=TRUE)
sumCOR$Grain.kg.ha[which(sumCOR$Grain.kg.ha=="NaN")]<-NA

ggplot(data=sumCOR, aes(x=sumCOR$Straw.kg.ha, y=sumCOR$Diff_IWG_LEG, color=sumCOR$Trt_code)) + 
  geom_point() +
  facet_grid(sumCOR$Location~sumCOR$fYear) +
  xlab("IWG Biomass (kg/ha)") + 
  ylab("15N Ratio Difference") + 
  ggtitle("IWG Biomass and IWG-Legume 15N Difference") +
  theme(panel.grid.minor=element_blank(),
        legend.title=element_blank())

cmod1 <- lm(Straw.kg.ha~Diff_IWG_LEG, data=LEG)
summary(cmod1)

cmod2<-lme(Straw.kg.ha~Diff_IWG_LEG, random=~1|Location/fYear, data=LEG, na.action=na.omit)
summary(cmod2)

ggplot(data=sumCOR, aes(x=sumCOR$Straw.kg.ha, y=sumCOR$Leg_Fraction, color=sumCOR$Trt_code)) + 
  geom_point() +
  facet_grid(sumCOR$Location~sumCOR$fYear) +
  xlab("IWG Biomass (kg/ha)") + 
  ylab("Legume-Dervied N Fraction") + 
  ggtitle("IWG Biomass and Legume N Fraction") +
  theme(panel.grid.minor=element_blank(),
        legend.title=element_blank())

cmod3 <- lm(Straw.kg.ha~Leg_Fraction, data=LEG)
summary(cmod3)

cmod4<-lme(Straw.kg.ha~Leg_Fraction, random=~1|Location/fYear, data=LEG, na.action=na.omit)
summary(cmod4)

cmod5 <- lm(Straw.kg.ha~Leg_Fraction, data=subset(LEG, Location=="Lamberton"))
summary(cmod5)

summary(lm(Straw.kg.ha~Leg_Fraction, data = SP17))
summary(lm(Straw.kg.ha~Leg_Fraction, data = SP18))
summary(lm(Straw.kg.ha~Leg_Fraction, data = SP19))
summary(lm(Straw.kg.ha~Leg_Fraction, data = Lam17))
summary(lm(Straw.kg.ha~Leg_Fraction, data = Lam18))
summary(lm(Straw.kg.ha~Leg_Fraction, data = Lam19))
summary(lm(Straw.kg.ha~Leg_Fraction, data = Ros17))
summary(lm(Straw.kg.ha~Leg_Fraction, data = Ros18))
summary(lm(Straw.kg.ha~Leg_Fraction, data = Ros19))



ggplot(data=sumCOR, aes(x=sumCOR$Leg.bio.kg.ha, y=sumCOR$Diff_IWG_LEG, color=sumCOR$Trt_code)) + 
  geom_point() +
  facet_grid(sumCOR$Location~sumCOR$fYear) +
  xlab("IWG Biomass (kg/ha)") + 
  ylab("15N Ratio Difference") + 
  ggtitle("Legume Biomass and IWG-Legume 15N Difference") +
  theme(panel.grid.minor=element_blank(),
        legend.title=element_blank())

cmod6 <- lm(Leg.bio.kg.ha~Diff_IWG_LEG, data=LEG)
summary(cmod6)

cmod7<-lme(Leg.bio.kg.ha~Diff_IWG_LEG, random=~1|Location/fYear, data=LEG, na.action=na.omit)
summary(cmod7)

cmod8 <- lm(Leg.bio.kg.ha~Diff_IWG_LEG, data=subset(LEG, Location=="Rosemount"))
summary(cmod8)

ggplot(data=sumCOR, aes(x=sumCOR$Leg.bio.kg.ha, y=sumCOR$Leg_Fraction, color=sumCOR$Trt_code)) + 
  geom_point() +
  facet_grid(sumCOR$Location~sumCOR$fYear) +
  xlab("IWG Biomass (kg/ha)") + 
  ylab("Legume-Dervied N Fraction") + 
  ggtitle("Legume Biomass and IWG Legume N Fraction") +
  theme(panel.grid.minor=element_blank(),
        legend.title=element_blank())

cmod9 <- lm(Leg.bio.kg.ha~Leg_Fraction, data=LEG)
summary(cmod9)

cmod10<-lme(Leg.bio.kg.ha~Leg_Fraction, random=~1|Location/fYear, data=LEG, na.action=na.omit)
summary(cmod10)

cmod11 <- lm(Leg.bio.kg.ha~Leg_Fraction, data=subset(LEG, Location=="Saint Paul"))
summary(cmod11)

#
cmod12 <- lm(IWG_CN~Straw.kg.ha, data=Rosdat)
summary(cmod12)

cmod13 <- lm(Straw.kg.ha~Leg.bio.kg.ha, data=SPdat)
summary(cmod13)
plot(Straw.kg.ha~Leg.bio.kg.ha, data=Lamdat)

#Grain and Legume Biomass
cmod14 <- lm(Grain.kg.ha~Leg.bio.kg.ha*Trt_code, data=LEG)
summary(cmod14)
plot(Grain.kg.ha~Leg.bio.kg.ha, data=LEG)

cmodA <- lm(Grain.kg.ha~Leg.bio.kg.ha*Trt_code, data=Rosdat)
summary(cmodA)
plot(Grain.kg.ha~Leg.bio.kg.ha, data=Rosdat, group=Trt_code)

cmodB <- lm(Grain.kg.ha~Leg.bio.kg.ha*Trt_code, data=Lamdat)
summary(cmodB)
plot(Grain.kg.ha~Leg.bio.kg.ha, data=Lamdat, group_by(Trt_code))

cmodC <- lm(Grain.kg.ha~Leg.bio.kg.ha*Trt_code, data=SPdat)
summary(cmodC)
plot(Grain.kg.ha~Leg.bio.kg.ha, data=SPdat)

ggplot(data=Rosdat, aes(x=Rosdat$Leg.bio.kg.ha, y=Rosdat$Grain.kg.ha, color=Rosdat$Trt_code)) + 
  geom_point() +
  facet_wrap(Rosdat$StandAge) +
  xlab("IWG Biomass (kg/ha)") + 
  ylab("Legume-Dervied N Fraction") + 
  ggtitle("Legume Biomass and IWG Legume N Fraction") +
  theme(panel.grid.minor=element_blank(),
        legend.title=element_blank())

ggplot(data=LEG, aes(x=LEG$Leg.bio.kg.ha, y=LEG$Grain.kg.ha, color=LEG$Trt_code)) + 
  geom_point() +
  facet_wrap(LEG$StandAge) +
  #facet_grid(sumCOR$Location~sumCOR$fYear) +
  xlab("IWG Biomass (kg/ha)") + 
  ylab("Legume-Dervied N Fraction") + 
  ggtitle("Legume Biomass and IWG Legume N Fraction") +
  theme(panel.grid.minor=element_blank(),
        legend.title=element_blank())

ggplot(data=LEG, aes(x=LEG$Leg.bio.kg.ha, y=LEG$Trt_code, color=LEG$Trt_code)) + 
  geom_point() +
  #facet_wrap(LEG$StandAge) +
  facet_grid(LEG$Location~LEG$StandAge) +
  xlab("Leg Biomass (kg/ha)") + 
  ylab("Treatment") + 
  ggtitle("Leg biomass by year") +
  theme(panel.grid.minor=element_blank(),
        legend.title=element_blank())


cmod15 <- lm(NTrans~Leg.bio.kg.ha, data=Lamdat)
summary(cmod15)

Legbio2 <- LEG$Leg.bio.kg.ha^2
quad.mod <- lm(Grain.kg.ha~Leg.bio.kg.ha + Legbio2, data=LEG)
summary(quad.mod)

plot(Grain.kg.ha~Leg.bio.kg.ha, data=LEG)

###Random Effects Models####
LEGStrawmod1<-lme(Straw.kg.ha ~ Trt_code * StandAge, data = LEG, random = ~1|fRep/Location, na.action = na.omit)
anova(LEGStrawmod1)
summary(LEGStrawmod1)
cld(emmeans(LEGStrawmod1, ~Trt_code))

LEGStrawmod2<-lme(Straw.kg.ha ~ Leg.bio.kg.ha*Trt_code, data = LEG, random = ~1|fRep/Location, na.action = na.omit)
anova(LEGStrawmod2)
summary(LEGStrawmod2)
cld(emmeans(LEGStrawmod2, ~Trt_code))

LEGGrainmod1<-lme(Grain.kg.ha ~ Trt_code, data = LEG, random = ~1|fRep/Location, na.action = na.omit)
anova(LEGGrainmod1)
summary(LEGGrainmod1)
cld(emmeans(LEGGrainmod1, ~Trt_code))

LEGLBiomod1<-lme(Leg.bio.kg.ha ~ Trt_code, data = LEG, random = ~1|fRep/Location, na.action = na.omit)
anova(LEGLBiomod1)
summary(LEGLBiomod1)
cld(emmeans(LEGLBiomod1, ~Trt_code))

LEGWBiomod1<-lme(Weed.bio.kg.ha ~ Trt_code, data = LEG, random = ~1|fRep/Location, na.action = na.omit)
anova(LEGWBiomod1)
summary(LEGWBiomod1)
cld(emmeans(LEGWBiomod1, ~Trt_code))

summarySE(LEG, "Leg.bio.kg.ha", c("Trt_code", "Location"), na.rm = TRUE)

#### Comparing intercrop grain yields to fertilized yields ####
sumgrain<-summarySE(LEG, "Grain.kg.ha", c("fYear", "Location","Trt_code"), na.rm=T)
env1<-subset(Rosdat, fStandAge=="3")
m1<-lme(Grain.kg.ha~Trt_code, data=subset(Rosdat, fStandAge=="3"), random = ~1|fRep, na.action = na.omit)
cld(emmeans(m1, ~Trt_code))
t.test(Grain.kg.ha~Trt_code, data=env1[which(env1$Trt_code==c("Control","Alfalfa")),])

summary(glht(m1, linfct=mcp(Trt_code = c("45 - Control = 0"))))
summary(glht(m1, mcp(Trt_code="Tukey")))
soildat<-soildat[which(soildat$trt!="22"),]
sumgrainL<-sumgrain[which(sumgrain$Trt_code!=c("Control","45","90")),]
sumgrainC<-sumgrain[which(sumgrain$Trt_code==c("Control")),]
sumgrainH<-sumgrain[which(sumgrain$Trt_code==c("90")),]
ggplot(sumgrainL, aes(y=Grain.kg.ha, x=Trt_code))+
  ggtitle("IWG grain yields from intercrops compared to Control")+
  facet_grid(Location~fYear)+
  geom_point()+
  geom_errorbar(aes(ymin=Grain.kg.ha-se, ymax=Grain.kg.ha+se), width=.2)+
  geom_hline(data=sumgrainC, aes(yintercept=Grain.kg.ha))+
  geom_hline(data=sumgrainC, aes(yintercept=Grain.kg.ha+se), linetype=2)+
  geom_hline(data=sumgrainC, aes(yintercept=Grain.kg.ha-se), linetype=2)

ggplot(sumgrainL, aes(y=Grain.kg.ha, x=Trt_code))+
  ggtitle("IWG grain yields from intercrops compared to 90 kg ha")+
  facet_grid(Location~fYear)+
  geom_point()+
  geom_errorbar(aes(ymin=Grain.kg.ha-se, ymax=Grain.kg.ha+se), width=.2)+
  geom_hline(data=sumgrainH, aes(yintercept=Grain.kg.ha))+
  geom_hline(data=sumgrainH, aes(yintercept=Grain.kg.ha+se), linetype=2)+
  geom_hline(data=sumgrainH, aes(yintercept=Grain.kg.ha-se), linetype=2)


