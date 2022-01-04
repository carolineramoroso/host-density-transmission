# re-approaching analysis
# using gams
library(ggplot2)
library(mgcv)
library(lmtest)
library(dplyr)
library(tidyr)
library(car)
library(tidymv)
library(plotrix)

theme_set(theme_bw())


#####***THROUGHOUT, OUTLINE OF PAPER IN BOLD AND WITH ASTERISKS***######

###############################################################################################
#####***POLLINATOR VISITS***#####

#PRIMARY VISITS ~ SPACING * DISEASE
poll_data = read.csv("pollinator visits_experiments 7 and 8.csv", header=T)
str(poll_data)

#######BEGIN DATA CLEANING STEPS#######
poll_data$Dis_or_H = factor(poll_data$Dis_or_H)
poll_data$o.Dis_or_H = ordered(poll_data$Dis_or_H)
poll_data$Rep.ID = factor(poll_data$Rep.ID)

poll_data$Spacing[poll_data$Spacing==0.30] = 0.33 #change spacing = 0.3 to 0.33
poll_data1 = poll_data[!poll_data$Spacing==0.67,] # remove 0.67 spacing

write.csv(poll_data1, "pollinator_visit_data_Rout.csv")  

#are there differences between replicates A and B? 
mod.rep = lm(log10(total_prim_visits_per_tube+0.042)~ Spacing + Rep.ID, data= poll_data1)
summary(mod.rep) #yes, and also sig diff if you ignore spacing

#analyze rep A and rep B separately
polldata_A = subset(poll_data1, Rep.ID==6)
sum(polldata_A$total_primvis)
polldata_B = subset(poll_data1, Rep.ID==7)
sum(polldata_B$total_primvis)

poll_data_means = poll_data1 %>% group_by(Spacing, Dis_or_H, o.Dis_or_H, Rep.ID) %>%
  summarize(Mean_vis_per_plant=mean(total_prim_visits_per_tube), mean_tubes=mean(Tubes))

polldatameans_A = subset(poll_data_means, Rep.ID==6)
polldatameans_B = subset(poll_data_means, Rep.ID==7)


#subset just the most abundant pollinator types for analysis - Rep B only 
poll_data_sub = polldata_B %>% select(c(Spacing, Dis_or_H, Butterflies_primvis_per_tube, Bees_primvis_per_tube,
                                        Hoverflies_primvis_per_tube, Houseflies_primvis_per_tube))
colnames(poll_data_sub)[3:6] = c("Butterflies", "Bees","Hoverflies","Houseflies")
poll_data_sub_gat = poll_data_sub %>% gather(Butterflies:Houseflies, key="Pollinator", value="Visits_per_tube")
poll_data_sub_gat$Pollinator = factor(poll_data_sub_gat$Pollinator, levels=c("Houseflies","Butterflies","Hoverflies","Bees"))

poll_data_sub_means = poll_data_sub_gat %>% group_by(Pollinator, Spacing) %>%
  summarize(Mean_vis_per_plant=mean(Visits_per_tube))


#######END DATA CLEANING STEPS#######

#######POLLINATOR VISITS ANALYSIS#######

#Total pollinator visits as a function of spacing and disease status

hist(log10(poll_data1$total_prim_visits_per_tube+0.042))


### *** TABLE S4 *** ### - use for stats
poll_data1$o.Rep.ID = ordered(poll_data1$Rep.ID)
poll_data1$o.Dis_or_H = ordered(poll_data1$Dis_or_H)

pol_gam1 = gam(log10(total_prim_visits_per_tube+0.042)~
                 s(Spacing, by=o.Rep.ID, k=5)  +
                 Rep.ID + Dis_or_H, data=poll_data1, method="REML")
summary(pol_gam1)
gam.check(pol_gam1)
anova(pol_gam1)


##### *** FIGURE 3 *** #####
plot_smooths(model = pol_gam1, series = Spacing, comparison = Dis_or_H, facet_terms=Rep.ID) + 
  stat_summary(fun.data = "mean_se", aes(x=Spacing, y=log10(total_prim_visits_per_tube+0.042), fill=Dis_or_H), 
               color="black", pch=21, alpha=0.7, data=poll_data1)+
  facet_wrap( ~factor(Rep.ID, levels=c(6,7),labels=c("Replicate A", "Replicate B")),  scales="free")+
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values=c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  ylab("log10(Total primary visits per plant+0.042)") + scale_x_reverse()


# analyzing replicates separately

# Replicate A
min(polldata_A$total_prim_visits_per_tube[polldata_A$total_prim_visits_per_tube!=0])
hist(log10(polldata_A$total_prim_visits_per_tube+0.042))

polA_gam1$o.Dis_or_H = ordered(polA_gam1$Dis_or_H)
polA_gam1 = gam(log10(total_prim_visits_per_tube+0.042)~s(Spacing, by=o.Dis_or_H, k=5) + Dis_or_H, data=polldata_A, method="REML")
summary(polA_gam1)
gam.check(polA_gam1)
anova(polA_gam1)

plot_smooths(model = polA_gam1, series = Spacing, comparison = Dis_or_H) + 
  stat_summary(fun.data = "mean_se", aes(x=Spacing, y=log10(total_prim_visits_per_tube+0.042), col=Dis_or_H, fill=Dis_or_H), 
               pch=21,alpha=0.7, data=polldata_A,position=position_dodge(width=0.05))+
  #geom_jitter(aes(x=Spacing, y=log10(total_prim_visits_per_tube+0.042), fill=Dis_or_H), 
  #           color="black", pch=21, alpha=0.7, size=2.5, data=polldata_A, width=0.01) + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))+
  ylab("Log10(Total primary visits per plant+0.042)") + xlab("Spacing (m)")+
  scale_x_reverse() + ggtitle("(a) Replicate A")

# Replicate B
min(polldata_B$total_prim_visits_per_tube[polldata_B$total_prim_visits_per_tube!=0])
hist(log10(polldata_B$total_prim_visits_per_tube+0.067))

polldata_B$o.Dis_or_H = ordered(polldata_B$Dis_or_H)
polB_gam1 = gam(log10(total_prim_visits_per_tube+0.067)~s(Spacing, by=o.Dis_or_H, k=5) + Dis_or_H, data=polldata_B, method="REML")
summary(polB_gam1)
gam.check(polB_gam1)
anova(polB_gam1)

plot_smooths(model = polB_gam1, series = Spacing, comparison = Dis_or_H) + 
  stat_summary(fun.data = "mean_se", aes(x=Spacing, y=log10(total_prim_visits_per_tube+0.042),col=Dis_or_H, fill=Dis_or_H), 
                pch=21,alpha=0.7, data=polldata_B, position=position_dodge(width=0.05))+
  #geom_jitter(aes(x=Spacing, y=log10(total_prim_visits_per_tube+0.067), fill=Dis_or_H), color="black", 
  #           pch=21, alpha=0.7, size=2.5, data=polldata_B, width=0.01) + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))+
  ylab("Log10(Total primary visits per plant+0.067)") + xlab("Spacing (m)")+
  scale_x_reverse() + ggtitle("(b) Replicate B")




#Visits by pollinator types - 4 most abundant, as a function of spacing
# only in rep B
poll_data_sub_gat
min(poll_data_sub_gat$Visits_per_tube[poll_data_sub_gat$Visits_per_tube!=0])
hist(log10(poll_data_sub_gat$Visits_per_tube+0.027))

### *** TABLE S5 *** ###
poll_data_sub_gat$Pollinator
polltypes_gam1 = gam(log10(Visits_per_tube+0.027)~s(Spacing, by=Pollinator, k=5) + Pollinator, data=poll_data_sub_gat, method="REML")
summary(polltypes_gam1)
gam.check(polltypes_gam1)
anova(polltypes_gam1)

### *** FIGURE 4 *** ### 
plot_smooths(model = polltypes_gam1, series = Spacing, comparison = Pollinator) + 
  stat_summary(fun.data="mean_se", 
               aes(x=Spacing, y=log10(Visits_per_tube+0.027), color=Pollinator), 
               alpha=0.7, data=poll_data_sub_gat, position= position_dodge(width=0.05))+
  ylab("log10(Primary visits per plant+0.027)") + scale_x_reverse() +
  scale_color_manual(values=c("#4daf4a","#984ea3","#e41a1c","#377eb8")) + 
  scale_fill_manual(values=c("#4daf4a","#984ea3","#e41a1c","#377eb8")) +
  scale_linetype_manual(values=c(1,1,1,1))


###############################################################################################
#####################################################################################
#####***FLORAL SPORE DEPOSITION***#####

# data from experiments 7 and 8
exp7 = read.csv("rep7-counts.csv", header=T)
exp8 = read.csv("rep8-counts.csv", header=T)

#######BEGIN DATA CLEANING STEPS#########
#drop a couple of columns and adjust column names 
#to make 7 and 8 datasets the same so they can be combined

new0.5 = read.csv("RepB_0.5mspacing.csv", header=T)
rep8 = read.csv("rep8-counts.csv", header=T)
rep7 = read.csv("rep7-counts.csv", header=T)

0.1456 # area per frame mm2
rep7$area_mm2 = rep7$frames * 0.1456
rep7$spores_mm2_new = rep7$total.count / rep7$area_mm2

rep8$area_mm2 = rep8$frames * 0.1456
rep8$spores_mm2_new = rep8$total.count / rep8$area_mm2

colnames(new0.5)[19] = "spores_mm2_new"
colnames(rep7)[23] = "column"

rep8no0.5 = rep8[rep8$spacing!=0.5,]

new0.5_sub = new0.5[,c("Sample","spacing","column","row","dist.from.dis","spores_mm2_new")]
rep8no0.5_sub = rep8no0.5[,c("Sample","spacing","column","row","dist.from.dis","spores_mm2_new")]
rep7_sub = rep7[,c("Sample","spacing","column","row","dist.from.dis","spores_mm2_new")]

compl_data8 = rbind(rep8no0.5_sub,new0.5_sub)
compl_data8$Rep = "B"
rep7_sub$Rep = "A"
compl_data = rbind(compl_data8, rep7_sub)
compl_data$Rep = factor(compl_data$Rep)

fl_sp_data = compl_data
###################
#account for the number of source plants
#fl_sp_data = read.csv("floral_spore_data_Rout1.csv", header=T)

fl_sp_data$nsources = NA
fl_sp_data$nsources[fl_sp_data$spacing==0.2] = 121 #0.2 spacing = 121 sources
fl_sp_data$nsources[fl_sp_data$spacing==0.3333] = 49 #0.33 = 49 
fl_sp_data$nsources[fl_sp_data$spacing==0.5] = 25 #0.5 = 25
fl_sp_data$nsources[fl_sp_data$spacing==0.6667] = 16 #0.67 = 16
fl_sp_data$nsources[fl_sp_data$spacing==1] = 9 #1 = 9 
fl_sp_data$nsources[fl_sp_data$spacing==2] = 4 #2 = 4

#get rid of some of the extra digits
fl_sp_data$spacing[fl_sp_data$spacing==0.3333] = 0.33 
fl_sp_data$spacing[fl_sp_data$spacing==0.6667] = 0.67 

#remove spacing = 0.67

fl_sp_data = fl_sp_data[!fl_sp_data$spacing==0.67,]

#get the total number of spores counted (not corrected per "frame" or per source)

#replace uncounted frames with mean of counted frames - 10 samples missing frames
#replaceNAs = rowSums(fl_sp_data[which(fl_sp_data$frames<9),4:12], na.rm=T) / fl_sp_data$frames[which(fl_sp_data$frames<9)] #get mean spore count of counted frames
#9 - fl_sp_data$frames[which(fl_sp_data$frames<9)] #figure out how many frames to replace in each row
#replaceNAslen = round(c(replaceNAs[1:8], rep(replaceNAs[9],times=4), rep(replaceNAs[10], times=3)),0) #make replacement list and round to nearest integer
#fl_sp_data[which(fl_sp_data$frames<9),4:12][is.na(fl_sp_data[which(fl_sp_data$frames<9),4:12])]=replaceNAslen
#fl_sp_data[which(fl_sp_data$frames<9),4:12] # no NAs

#fl_sp_data$count_tot = rowSums(fl_sp_data[,4:12], na.rm=T) 
#length(which(is.na(fl_sp_data$count_tot)))

fl_sp_data$Rep = factor(fl_sp_data$Rep)

write.csv(fl_sp_data, "floral_spore_data_Rout.csv")


# calculate mean spores deposited per target corrected for number of source plants
# i.e. for a whole plot, on average how many spores on each target per source?

#mean spores per source across spacing
floral_sporespersource = fl_sp_data %>% group_by(spacing, Rep, nsources) %>%
  summarize(mean_spores = mean(spores_mm2_new))
floral_sporespersource$mean_spores_per_source = floral_sporespersource$mean_spores/floral_sporespersource$nsources
floral_sporespersource

#get mean spores at each distance
#calculate std error of log transformed data (for plots)
floral_sp_meandist = fl_sp_data %>% group_by(spacing, Rep, dist.from.dis, nsources) %>%
  summarize(mean_spores = mean(spores_mm2_new), sterrlog = std.error(log10(spores_mm2_new+1)), n_samples = length(spores_mm2_new))

floral_sp_meandist$n_samples

plot(dist.from.dis~spacing, data = floral_sp_meandist)

#subset the more sampled spacing arrays 
#take the closest point, and then the closest points to 2m, 4m, and 6m 

fl_subset_mean= rbind(
  floral_sp_meandist[floral_sp_meandist$spacing==0.2 & (floral_sp_meandist$dist.from.dis==0.2|floral_sp_meandist$dist.from.dis==1.8|floral_sp_meandist$dist.from.dis==4.4|floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==0.33 & (floral_sp_meandist$dist.from.dis==0.3339|floral_sp_meandist$dist.from.dis==0.39|
                                                           floral_sp_meandist$dist.from.dis==2.0004|floral_sp_meandist$dist.from.dis==2.04|
                                                           floral_sp_meandist$dist.from.dis==4.0002|floral_sp_meandist$dist.from.dis==4.02|
                                                           floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==0.5 & (floral_sp_meandist$dist.from.dis==0.5 | floral_sp_meandist$dist.from.dis==2|floral_sp_meandist$dist.from.dis==4|floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==1 & (floral_sp_meandist$dist.from.dis==1|floral_sp_meandist$dist.from.dis==2|floral_sp_meandist$dist.from.dis==4|floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==2,]
)

plot(dist.from.dis~spacing, data = floral_sp_meandist)


#fl_sp_data_subset = rbind(
#  fl_sp_data[fl_sp_data$spacing==0.2 & (fl_sp_data$dist.from.dis==0.2|fl_sp_data$dist.from.dis==1.8|fl_sp_data$dist.from.dis==4.4|fl_sp_data$dist.from.dis==6),],
#  fl_sp_data[fl_sp_data$spacing==0.33 & (fl_sp_data$dist.from.dis==0.3339|fl_sp_data$dist.from.dis==0.39|
#                                                           fl_sp_data$dist.from.dis==2.0004|fl_sp_data$dist.from.dis==2.04|
#                                                           fl_sp_data$dist.from.dis==4.0002|fl_sp_data$dist.from.dis==4.02|
#                                                           fl_sp_data$dist.from.dis==6),],
#  fl_sp_data[fl_sp_data$spacing==0.5 & (fl_sp_data$dist.from.dis==2|fl_sp_data$dist.from.dis==4|fl_sp_data$dist.from.dis==6),],
#  fl_sp_data[fl_sp_data$spacing==1 & (fl_sp_data$dist.from.dis==1|fl_sp_data$dist.from.dis==2|fl_sp_data$dist.from.dis==4|fl_sp_data$dist.from.dis==6),],
#  fl_sp_data[fl_sp_data$spacing==2,])
#View(fl_sp_data_subset)

#fl_sp_data_subset_means = fl_sp_data_subset %>% group_by(spacing, Rep, dist.from.dis, nsources) %>%
#  summarize(mean_spores = mean(count_tot), n_samples = length(count_tot))
#fl_sp_data_subset_means

#######END OF DATA CLEANING STEPS#########

#######FLORAL SPORE ANALYSIS#######

#make spacing an ordered factor
fl_subset_mean$spacing_f = factor(fl_subset_mean$spacing, levels=c(2,1,0.5,0.33,0.2))
fl_subset_mean$o.spacing = ordered(fl_subset_mean$spacing_f, levels=c(2,1,0.5,0.33,0.2))
floral_sp_meandist$spacing_f = ordered(floral_sp_meandist$spacing, levels=c(2,1,0.5,0.33,0.2))
floral_sp_meandist$o.spacing = ordered(floral_sp_meandist$spacing_f, levels=c(2,1,0.5,0.33,0.2))

fl_subset_mean$log10mean_spores1=log10(fl_subset_mean$mean_spores+1)
floral_sp_meandist$log10mean_spores1=log10(floral_sp_meandist$mean_spores+1)

#ordered and regular factor for interpretation
### *** TABLE 1 *** ###
fl_gam1 = gam(log10mean_spores1~ s(dist.from.dis, by=o.spacing, k=5) + dist.from.dis + spacing_f + Rep, data=fl_subset_mean, weights = n_samples, method="REML")
summary(fl_gam1)
gam.check(fl_gam1)

anova(fl_gam1)


#
fl_gam1 = gam(log10mean_spores1~ s(dist.from.dis, by=o.spacing, k=5) + dist.from.dis + o.spacing + Rep, data=fl_subset_mean, weights = n_samples, method="REML")
summary(fl_gam1)

#just means - with shadow points
#plot_smooths(model = fl_gam1, series = dist.from.dis, comparison = Rep, facet_terms = o.spacing) + 
  #geom_point(data= fl_subset_mean, alpha=0.7, size=2.5, col="black",
  #           aes(x=dist.from.dis, y= log10(mean_spores+1), shape=Rep, fill=Rep)) +
  #geom_point(data= floral_sp_meandist, alpha=0.4, size=2.5, col="black",
  #           aes(x=dist.from.dis, y= log10(mean_spores+1), shape=Rep, fill=Rep)) +
  #facet_wrap(~factor(o.spacing,labels=c("Spacing=2m","Spacing=1m","Spacing=0.5m","Spacing=0.33m","Spacing=0.2m")),
  #           nrow=1)+
  #scale_shape_manual(values=c(21,23),name="Replicate", labels=c("A", "B")) +
  #scale_color_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  #scale_fill_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  #scale_linetype_manual(values=c("solid","dashed"), name="Replicate", labels=c("A", "B")) +
  #theme_bw() + ylab("Log10(Mean spores + 1)") + xlab("Distance from nearest diseased plant (m)")

######## *** FIGURE 5 *** ######

plot_smooths(model = fl_gam1, series = dist.from.dis, comparison = Rep, facet_terms = o.spacing) + 
  #all mean points in "shadow"
  geom_errorbar(data= floral_sp_meandist, alpha=0.6,
               aes(col=Rep, x=dist.from.dis, ymin= log10mean_spores1-sterrlog, ymax= log10mean_spores1+sterrlog)) +
  geom_point(data= floral_sp_meandist, alpha=0.6, size=2.5,
             aes(x=dist.from.dis, y= log10mean_spores1, shape=Rep,fill=Rep), col="gray40") +
  #subset of mean points that were analyzed
  geom_errorbar(data= fl_subset_mean, 
                aes(col=Rep, x=dist.from.dis, ymin= log10mean_spores1-sterrlog, ymax= log10mean_spores1+sterrlog)) +
  geom_point(data= fl_subset_mean, alpha=1, size=2.5,
             aes(x=dist.from.dis, y= log10mean_spores1, shape=Rep, fill=Rep), col="black") +
  facet_wrap(~factor(o.spacing,labels=c("Spacing=2m","Spacing=1m","Spacing=0.5m","Spacing=0.33m","Spacing=0.2m")),
             nrow=1)+
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("A", "B")) +
  scale_color_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_fill_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_linetype_manual(values=c("solid","dashed"), name="Replicate", labels=c("A", "B")) +
  theme_bw() + ylab("Log10(Mean spores + 1)") + xlab("Distance from nearest diseased plant (m)")




###############################################################################################
###############################################################################################
#####***AERIAL SPORE DEPOSITION***#####

### BEGIN DATA CLEANING STEPS ###

aerial8 = read.csv("Aerial data Expt 8_JA_1.2022.csv", header=T)
aerial7 = read.csv("rep 7 aerial.csv", header=T)

#1*clus1+10*clus10+100*clus100+1000*clus1000, then sum and log+1
<<<<<<< HEAD
#spore counts for rep 8 arrived at following this ^

aerial8$spores = aerial8$clus1 + 10*aerial8$clus10 + 100*aerial8$clus100 + 1000*aerial8$clus1000

=======
#spore counts for rep 8 arrived at following this 
>>>>>>> 9d2d99799a96b7195af049054366c2e7b93a7c48
aerial8_sub = aerial8[,c('plot', 'spacing', 'xx', 'type', 'spores')]
aerial8_sub$Experiment = c(8)

aerial7_sub = aerial7[,c('plot', 'spacing', 'xx', 'type', 'total_spores')]
colnames(aerial7_sub) = c('plot', 'spacing', 'xx', 'type', 'spores')
aerial7_sub$Experiment = c(7)
unique(aerial7_sub$spacing)

aerial = rbind(aerial7_sub, aerial8_sub)
aerial$spacing = round(aerial$spacing, 2)
aerial$xx = round(aerial$xx, 2)


aerial$nsources = NA

aerial$nsources[aerial$spacing==0.20] = 121 #0.2 spacing = 121 sources
aerial$nsources[aerial$spacing==0.33] = 49 #0.33 = 49 
aerial$nsources[aerial$spacing==0.50] = 25 #0.5 = 25
aerial$nsources[aerial$spacing== 1.00] = 9 #1 = 9 
aerial$nsources[aerial$spacing==2.00] = 4 #2 = 4
aerial$spores_per_source = aerial$spores / aerial$nsources
min(aerial$spores_per_source[aerial$spores_per_source != 0]) 

aerial$Experiment = factor(aerial$Experiment, levels=c(7,8), labels=c("A","B"))

<<<<<<< HEAD
write.csv(aerial, "aerial_spore_data_Rout.1.4.22.csv")
=======
write.csv(aerial, "aerial_spore_data_Rout.csv")
>>>>>>> 9d2d99799a96b7195af049054366c2e7b93a7c48

#"fold over" the aerial data to measure distance from disease in either direction
#take only the plates on the edge of the diseased -remove any D type that are not 6 or 8
plot(xx~spacing, aerial[aerial$type=="D",])
aerialnoD = aerial[aerial$type!="D",]
aerialDonly = aerial[aerial$type=="D",]
aerialDonly68 = aerialDonly[aerialDonly$xx==6 | aerialDonly$xx==8,]

aerial = rbind(aerialnoD, aerialDonly68)

aerial
aerial$dist.from.dis = NA
aerial$dist.from.dis[aerial$type=="D"] = 0
aerial$dist.from.dis[aerial$type=="H"] = 6-aerial$xx[aerial$type=="H"]
aerial$dist.from.dis[aerial$type=="E"] = aerial$xx[aerial$type=="E"]-8

aerial$dist.from.dis


plot(dist.from.dis~xx, aerial) #check
aerial$type = factor(aerial$type, levels=c("E","H","D"))
summary(lm( log10(spores+1)~Experiment + spacing + dist.from.dis + type, data=aerial))

#does it make a difference if the spores were near the healthy plants or in the empty space?
aerial.model = lm(log10(spores+1)~ dist.from.dis+type+spacing, data=aerial[aerial$type!="D",])
summary(aerial.model) # no... though if I leave in the ones within the diseased plants, D is different from H and E.

aerial$type.o = ordered(aerial$type)
gam_aer.eh = gam(log10(spores+1)~ s(xx, by=type.o) + type + spacing, data=aerial[aerial$type!="D",])
summary(gam_aer.eh) # no... though if I leave in the ones within the diseased plants, D is different from H and E.
anova(gam_aer.eh)


aerial2 = aerial

ggplot(data=aerial2, aes(x=xx, y=log10(spores+1))) + 
  stat_summary(aes(col=type), fun.data = "mean_se") + xlab("X position (m)") + ylab("Log10(Mean spore count+1)") +
  facet_wrap(~factor(spacing, levels=c(2,1,0.5,0.33,0.2), 
                     labels=c("Spacing=2m", "Spacing=1m","Spacing=0.5m","Spacing=0.33m","Spacing=0.2m")), nrow=1) +
  scale_color_manual(values=c("black","#78c679","#fa9fb5"), labels=c("Diseased","Empty","Healthy"), name="Surrounding Flowers")



#mean spores per source across spacing
aerial_sporespersource = aerial %>% group_by(spacing, Experiment, nsources) %>%
  summarize(mean_spores = mean(spores))
aerial_sporespersource$mean_spores_per_source = aerial_sporespersource$mean_spores/aerial_sporespersource$nsources
aerial_sporespersource


#mean spores at each distance

aerial_sp_meandist = aerial %>% group_by(spacing, Experiment, dist.from.dis, type, nsources) %>%
  summarize(mean_spores = mean(spores), sterrlog = std.error(log10(spores+1)), n_samples = length(spores))

aerial_sp_meandist$n_samples

plot(dist.from.dis~spacing, data = aerial_sp_meandist)

#subset the more sampled spacing arrays 
#take the closest point, and then the closest points to 2m, 4m, and 6m 
aer_subset_mean= aerial_sp_meandist[aerial_sp_meandist$dist.from.dis==0.00 |
                                      aerial_sp_meandist$dist.from.dis==1.00 | 
                                      aerial_sp_meandist$dist.from.dis==2.00 
                                    ,]
plot(dist.from.dis~spacing, data = aer_subset_mean)

summary(lm(mean_spores~dist.from.dis+type+spacing, data=aer_subset_mean))
summary(lm(mean_spores~type, data=subset(aer_subset_mean,type!="D"))) #no difference 



##### END DATA CLEANING STEPS#####
#########################################################################
#analysis
aer_subset_mean$spacing_f = factor(aer_subset_mean$spacing, levels=c(2,1,0.5,0.33,0.2))
aer_subset_mean$o.spacing = ordered(aer_subset_mean$spacing, levels=c(2,1,0.5,0.33,0.2))
#aer_subset_mean$o.spacing = ordered(aer_subset_mean$spacing, levels=c(0.2,0.33,0.5,1,2))

#log-transformed spore counts
aer_subset_mean$log10mean_spores1=log10(aer_subset_mean$mean_spores+1)
aerial_sp_meandist$log10mean_spores1=log10(aerial_sp_meandist$mean_spores+1)

#non-ordered factor for spacing for interpretation
aerialmeans_gam_log_wt = gam(log10mean_spores1~s(dist.from.dis, by= o.spacing, k=3) + dist.from.dis + Experiment + spacing_f, 
                             data=aer_subset_mean, method="REML", weights=n_samples)
gam.check(aerialmeans_gam_log_wt)
summary(aerialmeans_gam_log_wt)

###*** TABLE 2 ***###
anova(aerialmeans_gam_log_wt)

# ordered spacing factor for plotting 
aerialmeans_gam_log_wt = gam(log10mean_spores1~s(dist.from.dis, by= o.spacing, k=3) + dist.from.dis + Experiment + o.spacing, 
                             data=aer_subset_mean, method="REML", weights=n_samples)
summary(aerialmeans_gam_log_wt)

aerial_sp_meandist$o.spacing = ordered(aerial_sp_meandist$spacing, levels=c(2,1,0.5,0.33,0.2))

###### *** FIGURE 6 *** #####
# means with errorbars - with shadow points
plot_smooths(model = aerialmeans_gam_log_wt, series = dist.from.dis, comparison = Experiment, facet_terms = o.spacing) + 
  geom_errorbar(data= aer_subset_mean, width=0, alpha=1,
                aes(col=Experiment, x=dist.from.dis, ymin= log10mean_spores1-sterrlog, ymax= log10mean_spores1+sterrlog)) +
  geom_point(data= aer_subset_mean, col="black", alpha=1, size=2.5,
              aes(x=dist.from.dis, y= log10(mean_spores+1), 
                  shape=Experiment, fill=Experiment)) + 
  geom_errorbar(data= aerial_sp_meandist, width=0, alpha=0.6,
                aes(col=Experiment, x=dist.from.dis, ymin= log10mean_spores1-sterrlog, ymax= log10mean_spores1+sterrlog)) +
  geom_point(data= aerial_sp_meandist, alpha=0.6, size=2.5, col="black",
             aes(x=dist.from.dis, y= log10(mean_spores+1), 
                 shape=Experiment, fill=Experiment)) +
  facet_wrap(~ factor(o.spacing,labels=c("Spacing=2m","Spacing=1m","Spacing=0.5m","Spacing=0.33m","Spacing=0.2m")),
                     nrow=1) +
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("A", "B")) +
  scale_fill_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_color_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_linetype_manual(values=c("solid","dashed"), name="Replicate", labels=c("A", "B")) + 
  theme_bw() + ylab("Log10(Mean spores + 1)") + xlab("Distance from nearest diseased plant (m)")




#######################################
###############################################
#Aerial vs floral spore deposition against number of source plants

fl_subset_mean2 = fl_subset_mean[,c(1:5,7)]

#mean of spores on either side (grouping "type" E and H)
aer_subset_mean2 = aer_subset_mean %>% group_by(spacing, Experiment, dist.from.dis, nsources) %>%
  summarize(mean_spores = mean(mean_spores), n_samples=sum(n_samples))

colnames(aer_subset_mean2)[2] = "Rep" 
fl_subset_mean2$trans_mode = "Floral"
aer_subset_mean2$trans_mode = "Aerial"


both_types_subset = rbind(fl_subset_mean2,aer_subset_mean2)
tail(both_types_subset)

both_types_subset$Rep[both_types_subset$Rep=="Rep7"] = 7  
both_types_subset$Rep[both_types_subset$Rep=="Rep8"] = 8  
both_types_subset$Rep = factor(both_types_subset$Rep)
both_types_subset$trans_mode = factor(both_types_subset$trans_mode)
#both_types_subset$trans_mode = ordered(both_types_subset$trans_mode)


both_types_subset_2m = both_types_subset %>% filter(dist.from.dis==2 | dist.from.dis==2.0004 | dist.from.dis==2.04 | dist.from.dis==1.8)

ggplot(data=subset(both_types_subset_2m, Rep==7), aes(x=nsources, y=log10(mean_spores+1))) +
  geom_point(aes(col=trans_mode)) +  #facet_wrap(~Rep) + 
  geom_smooth(method="lm", formula=y~x, aes(col=trans_mode))


ggplot(data=subset(both_types_subset_2m, Rep==8), aes(x=nsources, y=log10(mean_spores+1))) +
  geom_point(aes(col=trans_mode)) +  #facet_wrap(~Rep) + 
  geom_smooth(method="lm", formula=y~x, aes(col=trans_mode))

both_types_subset_2m$o.trans_mode = ordered(both_types_subset_2m$trans_mode)

###*** TABLE S6 ***###
nsrc_gam =gam(log10(mean_spores+1) ~ s(nsources, by=o.trans_mode, k=5) + s(nsources, by=Rep, k=5) + trans_mode + Rep, data=both_types_subset_2m, weights = n_samples, method="REML")
summary(nsrc_gam)
anova(nsrc_gam)

nsrc_gam =gam(log10(mean_spores+1) ~ s(nsources, by=trans_mode, k=5) + s(nsources, by=Rep, k=5) + trans_mode + Rep, data=both_types_subset_2m, weights = n_samples, method="REML")


###### *** FIGURE 7 *** ######
plot_smooths(model = nsrc_gam, series = nsources, comparison = trans_mode, facet_terms = Rep) + 
  #mean-se of all of the subset of points
  #stat_summary(data= both_types_subset, aes(x=nsources, y= log10(mean_spores+1), fill=trans_mode, col=trans_mode), 
  #             fun.data= "mean_se", geom="pointrange", position=position_dodge(width=5), alpha=0.3) + 
  #just the 2 m points - which is what the analyses were done on... 
  geom_point(data= both_types_subset_2m, col="black", alpha=0.7, size=2.5, pch=21,
             aes(x=nsources, y= log10(mean_spores+1), fill=trans_mode)) +
  facet_wrap(~factor(Rep, labels=c("Replicate A","Replicate B")), scales="free") +
  xlab("Number of source diseased plants") + 
  ylab("Log10(Mean spore count + 1)") + 
  scale_color_manual(values=c("#4d9221","#c51b7d"), name="Transmission", labels=c("Aerial", "Floral")) +
  scale_fill_manual(values=c("#4d9221","#c51b7d"), name="Transmission", labels=c("Aerial", "Floral"))  +
  scale_linetype_manual(values=c("solid","dashed"), name="Transmission", labels=c("Aerial", "Floral")) 


###################################################################
### Control plots -- SPORES ###

# 
cont_spores = read.csv("control_spore counts_corrected.csv", header=T)
str(cont_spores)


cont_spores= cont_spores[,c(5:8,11:12)]

cont_spores$spacing[cont_spores$Plot=="A"] = 1
cont_spores$spacing[cont_spores$Plot=="B"] = 0.33
cont_spores$spacing[cont_spores$Plot=="C"] = 1
cont_spores$spacing[cont_spores$Plot=="D"] = 0.5
cont_spores$spacing[cont_spores$Plot=="E"] = 1
cont_spores$spacing[cont_spores$Plot=="F"] = 4 # single
cont_spores$spacing[cont_spores$Plot=="G"] = 4 # single

# 0.1685 mm2 per frame in 2019
cont_spores$count = as.numeric(cont_spores$count)
cont_spores$frames = as.numeric(cont_spores$frames)

cont_spores$area_mm2 = cont_spores$frames * 0.1685
cont_spores$spores_mm2_new = cont_spores$count /  cont_spores$area_mm2

#control spore deposition
cont_spores$x.pos*cont_spores$spacing

#"distance from disease" similar to experimental plots
cont_spores$dist.from.dis = cont_spores$x.pos*cont_spores$spacing
cont_spores$dist.from.dis[cont_spores$spacing==4] = 3
plot(spores_mm2_new~dist.from.dis, cont_spores)

cont_spores$dist.from.dis = round(cont_spores$dist.from.dis)

write.csv(cont_spores, "control_spore_data_Rout.csv")

#now get the subset of points and take means and SE at 0, 2, 4, 6 m 
cont_spores_meanpos = cont_spores %>% filter(dist.from.dis==0 | dist.from.dis==2 | dist.from.dis==4 | dist.from.dis==6 | spacing==4) %>%
  group_by(spacing, dist.from.dis) %>% summarize(mean_spores = mean(spores_total), sterrlog = std.error(log10(spores_total+1)), n_samples = length(spores_total))

subs_cont_spores_meanpos=cont_spores_meanpos
#subs_cont_spores_meanpos = subset(cont_spores_meanpos, spacing<4)
subs_cont_spores_meanpos$spacing_f = factor(subs_cont_spores_meanpos$spacing, levels=c(4,1,0.5,0.33))
subs_cont_spores_meanpos$o.spacing = ordered(subs_cont_spores_meanpos$spacing_f, levels=c(4,1,0.5,0.33))


cont_spores$spacing_f = factor(cont_spores$spacing, levels=c(4,1,0.5,0.33))
cont_spores$o.spacing = ordered(cont_spores$spacing_f, levels=c(4,1,0.5,0.33))


#just the means at a subset of locations (0, 2, 4, 6 m; or just the middle location for spacing=4m)
### ***TABLE S7*** ###
cont_gam1 = gam(log10(mean_spores+1)~ s(dist.from.dis, by=o.spacing, k=3)  + spacing_f, 
                data=subs_cont_spores_meanpos, weights = n_samples, method="REML")
summary(cont_gam1)
gam.check(cont_gam1)
anova(cont_gam1)

plot_smooths(model = cont_gam1, series = dist.from.dis, comparison = spacing_f)  
  
cont_spores_meanpos$log10mean_spores1 = log10(cont_spores_meanpos$mean_spores+1)

cont_gam1 = gam(log10(mean_spores+1)~ s(dist.from.dis, by=o.spacing, k=3)  + o.spacing, 
                data=subs_cont_spores_meanpos, weights = n_samples, method="REML")
###### ***Figure 8*** #####
ggplot(aes(y=log10(mean_spores+1), x=dist.from.dis), data=cont_spores_meanpos) + 
  facet_wrap(~factor(spacing, levels=c(4, 1, 0.5, 0.33),
                     labels=c("Spacing 4m+","Spacing 1m","Spacing 0.5m","Spacing 0.33m")), nrow=1) +
  geom_errorbar(col="gray20", aes(x=dist.from.dis, ymin= log10mean_spores1-sterrlog, ymax= log10mean_spores1+sterrlog), width=0) + 
  geom_point(col="gray20", fill="gray", pch=22, size=2.5) + 
  theme_bw() + xlab("Distance (m)") + ylab("Log10(Mean spore count+1)") +
  theme(legend.position = "none")


cont_spores_means = cont_spores %>% group_by(spacing) %>% summarize(mean_spores_perplant = mean(spores_total))
cont_spores_means$ntargets=NA
cont_spores_means$ntargets[cont_spores_means$spacing==0.33] = 133
cont_spores_means$ntargets[cont_spores_means$spacing==0.5] = 60
cont_spores_means$ntargets[cont_spores_means$spacing==1] = 18
cont_spores_means$ntargets[cont_spores_means$spacing==4] = 1


#################################################################################################

#### fan design

# data from the fan design experiment
fandata = read.csv("fan arrangement data CRA.csv", header=T)

#calculate per mm2 counts
fandata$count = as.numeric(fandata$count)
fandata$frames = as.numeric(fandata$frames)

fandata$area_mm2 = fandata$frames * 0.1685
fandata$spores_mm2_new = fandata$count /  fandata$area_mm2



library(spatstat.geom)

# get distance to nearest disease (min dist to disease)
#locations of diseased
ycoords = c(27.5004115, 17.1041181, 10.6380537, 6.6164292,	
            4.1151452, 2.55945, 1.591872, 0.9900785, 0.6157879) #just copied from excel
dis_coords = data.frame(xcor=rep(0,9), ycor=ycoords)

#get distances between sources and targets
distances_dis = crossdist.default(fandata$X_cor, fandata$Y_cor, dis_coords$xcor, dis_coords$ycor)

#distance to nearest diseased plant
fandata$min_dist_dis = apply(distances_dis, 1, FUN=min)
fandata$min_dist_dis = round(fandata$min_dist_dis, 1)
max(fandata$min_dist_dis)
fandata$count[fandata$min_dist_dis==8.9]

write.csv(fandata, "fan_data_Rout.csv")

fandata_means = fandata %>% group_by(min_dist_dis, Experiment_replicate) %>% summarize(mean_count= mean(spores_mm2_new), n_samples=length(spores_mm2_new))
#fandata_means$Experiment_replicate=factor(fandata_means$Experiment_replicate)
View(fandata_means)
fan_gam1 = gam(log10(spores_mm2_new+1)~ s(min_dist_dis, k=9), 
                data=fandata, method="REML")
summary(fan_gam1)
gam.check(fan_gam1)
anova(fan_gam1)

plot_smooths(model = fan_gam1, series = min_dist_dis) +
  #geom_point(data = fandata, aes(x=min_dist_dis, y=log10(count+1),col=Experiment_replicate)) +
  stat_summary(fun.data="mean_se", data = fandata, position=position_dodge(width=0.05),
               aes(x=min_dist_dis, y=log10(spores_mm2_new+1),fill=Experiment_replicate,col=Experiment_replicate, shape=Experiment_replicate))+
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("a", "b")) +
  scale_fill_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  scale_color_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  scale_linetype_manual(values=c("solid", "dashed"),name="Replicate", labels=c("a", "b")) +
  
  theme_bw() + ylab("Log10(Mean spore count+1)") + xlab("Distance from disease (m)") + xlim(0,10)


  #mean-se of all of the subset of points
  #stat_summary(data= both_types_subset, aes(x=nsources, y= log10(mean_spores+1), fill=trans_mode, col=trans_mode), 
  #             fun.data= "mean_se", geom="pointrange", position=position_dodge(width=5), alpha=0.3) + 
  #just the 2 m points - which is what the analyses were done on... 
  geom_point(data= both_types_subset_2m, col="black", alpha=0.7, size=2.5, pch=21,
             aes(x=nsources, y= log10(mean_spores+1), fill=trans_mode)) +
  facet_wrap(~factor(Rep, labels=c("Replicate A","Replicate B")), scales="free") +
  xlab("Number of source diseased plants") + 
  ylab("Log10(Mean spore count + 1)") + 
  scale_color_manual(values=c("#4d9221","#c51b7d"), name="Transmission", labels=c("Aerial", "Floral")) +
  scale_fill_manual(values=c("#4d9221","#c51b7d"), name="Transmission", labels=c("Aerial", "Floral"))  +
  scale_linetype_manual(values=c("solid","dashed"), name="Transmission", labels=c("Aerial", "Floral")) 



#Figure S7

ggplot(data=fandata_means, aes(x=min_dist_dis, y=mean_count)) +
  #geom_smooth(aes(y=prediction), method="lm", formula=y~x, col="gray40") +
  geom_point( aes(pch=Experiment_replicate, fill=Experiment_replicate), size=2.5, alpha=0.7) + theme_bw() + xlab("Minimum distance from disease (m)") +
  ylab("Mean spore count") +
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("a", "b")) +
  scale_fill_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  scale_color_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  theme_bw() + ylab("Log10(Mean spore count+1)") + xlab("Distance from disease (m)")+ xlim(0,10)


max(fandata_means$min_dist_dis)
(fandata$count[fandata$min_dist_dis==8.93])

ggplot(data=fandata_means, aes(x=min_dist_dis, y=mean_count)) +
  #geom_smooth(aes(y=prediction), method="lm", formula=y~x, col="gray40") +
  geom_point( data = fandata, position=position_dodge(width=0.05), alpha=0.7,size=2.5,
               aes(x=min_dist_dis, y=count,fill=Experiment_replicate,col=Experiment_replicate, shape=Experiment_replicate))+
  ylab("Mean spore count") +
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("a", "b")) +
  scale_fill_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  scale_color_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  theme_bw() + ylab("Log10(Mean spore count+1)") + xlab("Distance from disease (m)")+ xlim(0,10)





#
ggplot(data=fandata, aes(x=min_dist_dis, y=log10(count+1))) +
  geom_jitter(fill="#78c679", pch=21, size=2.5, alpha=0.7, width=0.2) + theme_bw() + xlab("Minimum distance to disease") +
  ylab("Log10(Spore count+1)")+ xlim(0,10)







################# Secondary visits by pollinators 


########## pollinator secondary visits
poll_data = read.csv("pollinator visits_experiments 7 and 8.csv", header=T)

poll_data$Spacing[poll_data$Spacing == 0.30]= 0.33
poll_data = poll_data[!poll_data$Spacing==0.67,] # remove 0.67 spacing

poll_data$total_secvis = NA

str(poll_data)
for(i in 1:length(poll_data$total_secvis)){
  poll_data$total_secvis[i] = length(which(poll_data[i,11:59]>1))
}

sum(poll_data$total_primvis[poll_data$Rep.ID==6])
sum(poll_data$total_primvis[poll_data$Rep.ID==7])
sum(poll_data$total_primvis)
sum(poll_data$total_secvis)
sum(poll_data$total_secvis)/sum(poll_data$total_primvis)

poll_data_secvis = poll_data[,c(1:10,75:77)]

poll_data_secvis$prob_sec_given_prim = poll_data_secvis$total_secvis / poll_data_secvis$total_primvis
poll_data_secvis$prob_sec_perobsplant = poll_data_secvis$prob_sec_given_prim / poll_data_secvis$Tubes

write.csv(poll_data_secvis, "pollinator_second_visits_Rout.csv")


vis_primsec_mean = poll_data_secvis %>% group_by(Spacing, Rep.ID, Dis_or_H) %>% 
  summarize(mean_probsecvis= mean(prob_sec_given_prim, na.rm=T), mean_probsecvispp = mean(prob_sec_perobsplant, na.rm=T))
plot(mean_probsecvispp~Spacing, data= vis_primsec_mean)
plot(mean_probsecvis~Spacing, data= vis_primsec_mean)


secvis_gam1 = gam(mean_probsecvispp~s(Spacing, by=Dis_or_H, k=5) +
                   + Dis_or_H + Rep.ID, data=vis_primsec_mean, method="REML")
gam.check(secvis_gam1)
summary(secvis_gam1)
anova(secvis_gam1)

plot_smooths(model = secvis_gam1, series=Spacing, comparison=Dis_or_H, facet_terms=Rep.ID) + 
  geom_point(data=vis_primsec_mean, aes(x=Spacing, y=mean_probsecvispp, col=Dis_or_H)) + 
  facet_wrap(~factor(Rep.ID, labels =c("Replicate A","Replicate B")), scales="free") + 
  scale_x_reverse() + ylab("Probability of secondary visit given primary") + xlab("Spacing (m)") + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))


secvis_gam2 = gam(mean_probsecvis~s(Spacing, by=Dis_or_H, k=5) +
                    + Dis_or_H + Rep.ID, data=vis_primsec_mean, method="REML")
gam.check(secvis_gam2)
summary(secvis_gam2)
anova(secvis_gam2) #

plot_smooths(model = secvis_gam2, series=Spacing, comparison=Dis_or_H, facet_terms=Rep.ID) + 
  geom_point(data=vis_primsec_mean, aes(x=Spacing, y=mean_probsecvis, col=Dis_or_H)) + 
  facet_wrap(~factor(Rep.ID, labels =c("Replicate A","Replicate B")), scales="free") + 
  scale_x_reverse() + ylab("Probability of secondary visit given primary") + xlab("Spacing (m)") + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))



sec_vis_gat = gather(second_visits, key=PollType, value=VisPerTube, Butterflies:Total)
sec_vis_gat

sum(sec_vis_gat$VisPerTube[-which(sec_vis_gat$PollType=="Total")] * sec_vis_gat$Tubes[-which(sec_vis_gat$PollType=="Total")])
sum(sec_vis_gat$VisPerTube[which(sec_vis_gat$PollType=="Total")] * sec_vis_gat$Tubes[which(sec_vis_gat$PollType=="Total")])

sec_vis_gat$PollType = factor(sec_vis_gat$PollType, levels=c("Beetles", "Bees", "Butterflies", "Medium flies", "Hoverflies","Small flies", "Other", "Total"))

sec_vis_mean = sec_vis_gat %>% filter(PollType == "Total", Spacing != 0.67) %>%
  group_by(Spacing, Rep.ID, Dis_or_H) %>% summarize(mean_secvisperplant = mean(VisPerTube, na.rm=T))

poll_data_means
sec_vis_mean

sec_vis_meansperprim=merge(poll_data_means, sec_vis_mean, by=c("Rep.ID", "Dis_or_H", "Spacing"), all=T)
sec_vis_meansperprim$sec_per_prim = sec_vis_meansperprim$mean_secvisperplant / sec_vis_meansperprim$Mean_vis_per_plant
sec_vis_meansperprim$sec_per_prim_pertub = sec_vis_meansperprim$sec_per_prim /sec_vis_meansperprim$mean_tubes

sec_vis_meansperprim$o.Dis_or_H = ordered(sec_vis_meansperprim$Dis_or_H)
sec_vis_meansperprim$o.Rep.ID = ordered(sec_vis_meansperprim$Rep.ID)

#start here
secvis_gam3 = gam(sec_per_prim~s(Spacing, by=o.Dis_or_H, k=5) +s(Spacing, by=o.Rep.ID, k=5) +
                    + Dis_or_H + Rep.ID, data=sec_vis_meansperprim, method="REML")

secvis_gam3 = gam(sec_per_prim~s(Spacing, by=Dis_or_H, k=5) +s(Spacing, by=Rep.ID, k=5) +
                 + Dis_or_H + Rep.ID, data=sec_vis_meansperprim, method="REML")
gam.check(secvis_gam3)
summary(secvis_gam3)
anova(secvis_gam3)

#figure S2-A
plot_smooths(model = secvis_gam3, series=Spacing, comparison=Dis_or_H, facet_terms=Rep.ID) + 
  geom_point(data=sec_vis_meansperprim, aes(x=Spacing, y=sec_per_prim, col=Dis_or_H)) + 
  facet_wrap(~factor(Rep.ID, labels =c("Replicate A","Replicate B")), scales="free") + 
  scale_x_reverse() + ylab("Probability of secondary visit given primary") + xlab("Spacing (m)") + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))



secvis_gam4 = gam(sec_per_prim_pertub~s(Spacing, by=o.Dis_or_H, k=5) +s(Spacing, by=o.Rep.ID, k=5) +
                    + Dis_or_H + Rep.ID, data=sec_vis_meansperprim, method="REML")

secvis_gam4 = gam(sec_per_prim_pertub~s(Spacing, by=Dis_or_H, k=5) +s(Spacing, by=Rep.ID, k=5) +
                   + Dis_or_H + Rep.ID, data=sec_vis_meansperprim, method="REML")
gam.check(secvis_gam4)
summary(secvis_gam4)
anova(secvis_gam4)

#figure S2-B
plot_smooths(model = secvis_gam4, series=Spacing, comparison=Dis_or_H, facet_terms=Rep.ID) + 
  geom_point(data=sec_vis_meansperprim, aes(x=Spacing, y=sec_per_prim_pertub, col=Dis_or_H)) + 
  facet_wrap(~factor(Rep.ID, labels =c("Replicate A","Replicate B")), scales="free") + 
  scale_x_reverse() + ylab("Probability of secondary visit given primary - per plant") + xlab("Spacing (m)") + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))
  


secvis_gam = gam(sec_per_prim~s(Spacing, by=Dis_or_H, k=5) 
                 + Dis_or_H + Rep.ID, data=sec_vis_meansperprim, method="REML")
gam.check(secvis_gam)
summary(secvis_gam)



########

#compare the two replicates and control
cont_spores_meansp = cont_spores %>% group_by(spacing) %>% summarize(Rep="control", mean_spores_spac = mean(spores_mm2_new))

expt_spores_meansp = fl_sp_data %>% group_by(spacing, Rep) %>% summarize(mean_spores_spac=mean(spores_mm2_new))

meanspores_mm = rbind(cont_spores_meansp,expt_spores_meansp)
colnames(meanspores_mm)=c("Spacing","Replicate","mean_spores_spac")

meanpollvis = poll_data %>% subset(Disease_status=="H")%>% group_by(Spacing, Replicate) %>% summarize(meanvispertube=mean(Total_per_tube))

contpollvis = read.csv("control array pollinator visitations.csv", header=T) %>% group_by(Plot, date.time, Density) %>%
  summarize(Total_Polls=sum(Pollinators), Tubes=mean(Num.plants))
contpollvis$Vis_per_tube = contpollvis$Total_Polls / contpollvis$Tubes

contpollvis2 = contpollvis %>% group_by(Density) %>% summarize(Replicate="control",meanvispertube = mean(Vis_per_tube)) 
colnames(contpollvis2)[1] = "Spacing"
contpollvis2$Spacing[contpollvis2$Spacing=="SINGLE"] = 4
contpollvis2$Spacing = as.numeric(contpollvis2$Spacing)

polvis_comb = rbind(meanpollvis, contpollvis2)

pol_spore_comb = merge(polvis_comb,meanspores_mm, by=c("Spacing","Replicate"))


ggplot(data=pol_spore_comb, aes(x=meanvispertube, y=mean_spores_spac, col=Replicate)) +
  geom_point() 




