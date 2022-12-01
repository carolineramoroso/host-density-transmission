#test 123
library(ggplot2)
library(mgcv)
library(lmtest)
library(dplyr)
library(tidyr)
library(car)
library(tidymv)
library(plotrix)

theme_set(theme_bw())
setwd("~/Git/Host density shapes transmission mode/host denisty data")


#####***POLLINATOR VISITATION***#####
poll_data = read.csv("pollinator_primary_visits.csv", header=T)

poll_data %>% group_by(Replicate) %>% summarize(tot_prim_vis = sum(Total), tot_tubes = sum(Tubes), tot_hrs = sum(Obs_duration)/60)

poll_data$Disease_status = factor(poll_data$Disease_status)
poll_data$Replicate = factor(poll_data$Replicate)

#### Overall primary visitation rate as a function of density ####

#are there differences between replicates A and B? 
mod.rep = lm(log10(Total_per_tube+0.042)~ Spacing + Replicate, data= poll_data)
summary(mod.rep) #yes, and also sig diff if you ignore spacing

#### TABLE S4 ####
#Generalized additive model of primary pollinator visitations, 
#where the visitation rate is the number of pollinators of all 
#types observed per plant during each observation period.

poll_data$o.Replicate = ordered(poll_data$Replicate)
poll_data$o.Disease_status = ordered(poll_data$Disease_status)

pol_gam1 = gam(log10(Total_per_tube+0.042) ~ s(Spacing, by=o.Replicate, k=5) + Disease_status + Replicate, 
               data=poll_data, method="REML")
summary(pol_gam1)
gam.check(pol_gam1)
anova(pol_gam1)

#### FIGURE S4 ####
#Rate of primary visitation by all pollinator types per observed 
#floral tube containing two Dianthus pavonius flowers that were either 
#healthy or diseased and arranged at different densities.

#use unordered factor for plotting
pol_gam1 = gam(log10(Total_per_tube+0.042) ~ s(Spacing, by=Replicate, k=5) + Disease_status + Replicate, 
               data=poll_data, method="REML")

plot_smooths(model = pol_gam1, series = Spacing, comparison = Disease_status, facet_terms=Replicate) + 
  stat_summary(fun.data = "mean_se", aes(x=Spacing, y=log10(Total_per_tube+0.042), fill=Disease_status), 
               color="black", pch=21, alpha=0.7, data=poll_data)+
  facet_wrap( ~factor(Replicate, labels=c("Replicate A", "Replicate B")),  scales="free")+
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values=c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  ylab("log10(Total primary visits per plant+0.042)") + scale_x_reverse()


#### Pollinator Types ####

#per-tube measurements for each poll type
poll_data$Butterflies_pt = poll_data$Butterflies/poll_data$Tubes
poll_data$Bees_pt = poll_data$Bees/poll_data$Tubes
poll_data$Hoverflies_pt = poll_data$Hoverflies/poll_data$Tubes
poll_data$Houseflies_pt = poll_data$Houseflies/poll_data$Tubes

#subset just the most abundant pollinator types for analysis - Rep B only 
poll_data_sub = poll_data %>% subset(Replicate == "B") %>%
  select(c(Spacing, Disease_status, Butterflies_pt, Bees_pt, Hoverflies_pt, Houseflies_pt))
colnames(poll_data_sub)[3:6] = c("Butterflies", "Bees","Hoverflies","Flies")

#rearrange columns to rows
poll_data_sub_gat = poll_data_sub %>% gather(Butterflies:Flies, key="Pollinator", value="Visits_per_tube")
poll_data_sub_gat$Pollinator = factor(poll_data_sub_gat$Pollinator, levels=c("Flies","Butterflies","Hoverflies","Bees"))

#### TABLE S5 ####
#Results from a generalized additive model of primary pollinator 
#visitations by the four  most common pollinators in experimental 
#replicate B, which had a large enough sample size of pollinators 
#to examine patterns across types. 

min(poll_data_sub_gat$Visits_per_tube[poll_data_sub_gat$Visits_per_tube!=0])
polltypes_gam1 = gam(log10(Visits_per_tube+0.027)~s(Spacing, by=Pollinator, k=5) + Pollinator, data=poll_data_sub_gat, method="REML")
summary(polltypes_gam1)
gam.check(polltypes_gam1)
anova(polltypes_gam1)

#### FIGURE 2 #### 
#Mean and standard error of visitation rates of the four most common 
#pollinator types in Replicate B across plant densities

plot_smooths(model = polltypes_gam1, series = Spacing, comparison = Pollinator) + 
  stat_summary(fun.data="mean_se", 
               aes(x=Spacing, y=log10(Visits_per_tube+0.027), color=Pollinator), 
               alpha=0.7, data=poll_data_sub_gat, position= position_dodge(width=0.05))+
  ylab("Log10(Primary visits per plant+0.027)") + scale_x_reverse() +
  xlab("Spacing (m)") +
  scale_color_manual(values=c("#4daf4a","#984ea3","#e41a1c","#377eb8")) + 
  scale_fill_manual(values=c("#4daf4a","#984ea3","#e41a1c","#377eb8")) +
  scale_linetype_manual(values=c(1,1,1,1)) + theme_bw()


#### Secondary pollinator visitations ####

second_visits = read.csv("pollinator_secondary_visits.csv", header=T)

sum(second_visits$Primary)
sum(second_visits$Secondary)
sum(second_visits$Secondary) / sum(second_visits$Primary)

second_vis_means = second_visits %>% group_by(Spacing, Replicate, Disease_status) %>% 
  summarize(mean_probsecvis= mean(Second_per_prim, na.rm=T), mean_probsecvispp = mean(Second_per_prim_per_plant, na.rm=T))

second_vis_means$Disease_status = factor(second_vis_means$Disease_status)
second_vis_means$o.Disease_status = ordered(second_vis_means$Disease_status)
second_vis_means$Replicate = factor(second_vis_means$Replicate)
second_vis_means$o.Replicate = ordered(second_vis_means$Replicate)

#### TABLE S6 - A ####
# Generalized additive model of secondary pollinator visitations. 
#A. The secondary visitation rate is the probability of a second visit given a primary visit, 
#calculated as the mean frequency of secondary visits per observed plant at each spacing divided 
#by the frequency of primary visits at each spacing.
secvis_gam1 = gam(mean_probsecvis~s(Spacing, by=o.Disease_status, k=5) +s(Spacing, by=o.Replicate, k=5) +
                    + Disease_status + Replicate, data=second_vis_means, method="REML")

gam.check(secvis_gam1)
summary(secvis_gam1)
anova(secvis_gam1)

#### FIGURE S5 - a ####
#Secondary pollinator visitations. Probability of secondary visit given a primary visit was 
#calculated by dividing the secondary visitation rate per observed plant by the 
#mean primary visitation rate per observed plant.

#remove ordered factors for purpose of plotting
secvis_gam1 = gam(mean_probsecvis~s(Spacing, by=Disease_status, k=5) +s(Spacing, by=Replicate, k=5) +
                    + Disease_status + Replicate, data=second_vis_means, method="REML")

plot_smooths(model = secvis_gam1, series=Spacing, comparison=Disease_status, facet_terms=Replicate) + 
  geom_point(data=second_vis_means, aes(x=Spacing, y=mean_probsecvis, col=Disease_status)) + 
  facet_wrap(~factor(Replicate, labels =c("Replicate A","Replicate B")), scales="free") + 
  scale_x_reverse() + ylab("Probability of secondary visit given primary") + xlab("Spacing (m)") + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))


#### TABLE S6 - B ####
#Results from generalized additive models of secondary pollinator visitations. 
#Corrected for the mean number of observed plants at each spacing.
secvis_gam2 = gam(mean_probsecvispp ~ s(Spacing, by=o.Disease_status, k=5) + s(Spacing, by=o.Replicate, k=5) +
                    + Disease_status + Replicate, data=second_vis_means, method="REML")

gam.check(secvis_gam2)
summary(secvis_gam2)
anova(secvis_gam2)


#### FIGURE S5 - b ####
#Secondary pollinator visitations. 
#corrected for the differences in the overall number of plants observed by 
#dividing this probability by the mean number of plants observed at each density

#remove ordered factors for purpose of plotting
secvis_gam2 = gam(mean_probsecvispp~s(Spacing, by=Disease_status, k=5) +s(Spacing, by=Replicate, k=5) +
                    + Disease_status + Replicate, data=second_vis_means, method="REML")

plot_smooths(model = secvis_gam2, series=Spacing, comparison=Disease_status, facet_terms=Replicate) + 
  geom_point(data=second_vis_means, aes(x=Spacing, y=mean_probsecvispp, col=Disease_status)) + 
  facet_wrap(~factor(Replicate, labels =c("Replicate A","Replicate B")), scales="free") + 
  scale_x_reverse() + ylab("Probability of secondary visit given primary - per plant") + xlab("Spacing (m)") + 
  scale_fill_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) +
  scale_color_manual(values=c("gray50","#c51b7d"), name="Plant status", labels=c("Diseased", "Healthy")) + 
  scale_linetype_manual(values =c("dashed","solid"), name="Plant status", labels=c("Diseased", "Healthy"))





#####***FLORAL SPORE DEPOSITION***#####

# data from experiments 7 and 8
fl_sp_data=read.csv("floral_spore_data.csv", header=T)
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

#first use the entire dataset

#make spacing an ordered factor
floral_sp_meandist$spacing_f = ordered(floral_sp_meandist$spacing, levels=c(2,1,0.5,0.33,0.2))
floral_sp_meandist$o.spacing = ordered(floral_sp_meandist$spacing_f, levels=c(2,1,0.5,0.33,0.2))

floral_sp_meandist$log10mean_spores1=log10(floral_sp_meandist$mean_spores+1)

#### TABLE 1A #### **************
#generalized additive model indicating that the relationship between 
#distance from the source of spores (in meters) and the number of spores deposited 
#on flowers is mediated by plant density.

fl_gam_entire = gam(log10mean_spores1~ s(dist.from.dis, by=o.spacing, k=5) + dist.from.dis + spacing_f + Rep, data=floral_sp_meandist, weights = n_samples, method="REML")
summary(fl_gam_entire)
gam.check(fl_gam_entire)
anova(fl_gam_entire)

#### FIGURE 3 - a ####
#Mean and standard error of the number of floral spores deposited on target plants by 
#vectors (pollinators) (a) 

fl_gam_entire = gam(log10mean_spores1~ s(dist.from.dis, by=o.spacing, k=5) + dist.from.dis + o.spacing + Rep, data=floral_sp_meandist, weights = n_samples, method="REML")
summary(fl_gam_entire)

plot_smooths(model = fl_gam_entire, series = dist.from.dis, comparison = Rep, facet_terms = spacing) + 
  geom_errorbar(data= floral_sp_meandist,
                aes(col=Rep, x=dist.from.dis, ymin= log10mean_spores1-sterrlog, ymax= log10mean_spores1+sterrlog)) +
  geom_point(data= floral_sp_meandist, size=2.5,
             aes(x=dist.from.dis, y= log10mean_spores1, shape=Rep,fill=Rep), col="gray40") +
  facet_wrap(~factor(o.spacing,labels=c("Spacing=2m","Spacing=1m","Spacing=0.5m","Spacing=0.33m","Spacing=0.2m")),
             nrow=1)+
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("A", "B")) +
  scale_color_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_fill_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_linetype_manual(values=c("solid","dashed"), name="Replicate", labels=c("A", "B")) +
  theme_bw() + ylab("Log10(Mean spores + 1)") + xlab("Distance from nearest diseased plant (m)")


### now consider that the denser plots have larger sample sizes, and it might just be 
### easier to pick up a non-linear shape when there is a higher density of points sampled along the
### x-axis. to adjust:
#subset the more sampled spacing arrays 
#take the closest point, and then the closest points to 2m, 4m, and 6m 
# (note, the results do not change, we are just being more conservative)

fl_subset_mean= rbind(
  floral_sp_meandist[floral_sp_meandist$spacing==0.2 & (floral_sp_meandist$dist.from.dis==0.2|floral_sp_meandist$dist.from.dis==1.8|floral_sp_meandist$dist.from.dis==4.4|floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==0.33 & (floral_sp_meandist$dist.from.dis==0.3339|floral_sp_meandist$dist.from.dis==0.39|
                                                           floral_sp_meandist$dist.from.dis==2.0004|floral_sp_meandist$dist.from.dis==2.04|
                                                           floral_sp_meandist$dist.from.dis==4.0002|floral_sp_meandist$dist.from.dis==4.02|
                                                           floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==0.5 & (floral_sp_meandist$dist.from.dis==0.5|floral_sp_meandist$dist.from.dis==2|floral_sp_meandist$dist.from.dis==4|floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==1 & (floral_sp_meandist$dist.from.dis==1|floral_sp_meandist$dist.from.dis==2|floral_sp_meandist$dist.from.dis==4|floral_sp_meandist$dist.from.dis==6),],
  floral_sp_meandist[floral_sp_meandist$spacing==2,]
)


#make spacing an ordered factor
fl_subset_mean$spacing_f = factor(fl_subset_mean$spacing, levels=c(2,1,0.5,0.33,0.2))
fl_subset_mean$o.spacing = ordered(fl_subset_mean$spacing_f, levels=c(2,1,0.5,0.33,0.2))
floral_sp_meandist$spacing_f = ordered(floral_sp_meandist$spacing, levels=c(2,1,0.5,0.33,0.2))
floral_sp_meandist$o.spacing = ordered(floral_sp_meandist$spacing_f, levels=c(2,1,0.5,0.33,0.2))

fl_subset_mean$log10mean_spores1=log10(fl_subset_mean$mean_spores+1)
floral_sp_meandist$log10mean_spores1=log10(floral_sp_meandist$mean_spores+1)


#### TABLE S7 ####
#Results from a generalized additive model indicating that the relationship between 
#distance from the source of spores (in meters) and the number of spores deposited 
#on flowers is mediated by plant density.

fl_gam1 = gam(log10mean_spores1~ s(dist.from.dis, by=o.spacing, k=5) + dist.from.dis + spacing_f + Rep, data=fl_subset_mean, weights = n_samples, method="REML")
summary(fl_gam1)
gam.check(fl_gam1)
anova(fl_gam1)



#####***AERIAL SPORE DEPOSITION***#####
aerial = read.csv("aerial_spore_data_4.6.22.csv", header=T)
#"fold over" the aerial data to measure distance from disease in either direction
#take only the plates on the edge of the diseased -remove any D type that are not 6 or 8

aerial
aerial$dist.from.dis = NA
aerial$dist.from.dis[aerial$type=="D"] = 0
aerial$dist.from.dis[aerial$type=="H"] = 6-aerial$xx[aerial$type=="H"]
aerial$dist.from.dis[aerial$type=="E"] = aerial$xx[aerial$type=="E"]-8

aerial$type = factor(aerial$type, levels=c("D","E","H"))

aerialA = subset(aerial, Experiment=="A")
aerialB = subset(aerial, Experiment=="B")

#### TABLE S8 ####
#chi square test 

repA_clumps = aerialA %>% group_by(type) %>% 
  summarize(Clumps10 = sum(clus1,clus10), Clumps100.1000 = sum(clus100,clus1000))
repA_clumps

A.D = repA_clumps[1,2:3] #clumps under diseased plants
A.E = repA_clumps[2,2:3] #clumps under no plants
A.H = repA_clumps[3,2:3] #clumps under healthy plants

A.nonD = A.H + A.E
nums = paste(c(A.D,A.nonD))
nums = as.numeric(nums)
dat.A= matrix(nums, nrow=2, ncol=2)
dat.A
chisq.test(x=dat.A)


repB_clumps = aerialB %>% group_by(type) %>% 
  summarize(Clumps10 = sum(clus1,clus10), Clumps100.1000 = sum(clus100,clus1000))
repB_clumps

B.D = repB_clumps[1,2:3] #clumps under diseased plants
B.E = repB_clumps[2,2:3] #clumps under no plants
B.H = repB_clumps[3,2:3] #clumps under healthy plants

B.nonD = B.H + B.E
nums = paste(c(B.D,B.nonD))
nums = as.numeric(nums)
dat.B= matrix(nums, nrow=2, ncol=2)
dat.B
chisq.test(x=dat.B)


###########
#### TABLE S9 ####
#no difference detected between aerial spore traps that were placed among healthy 
#flowers and in empty space.  
aerialA_noD = subset(aerialA, type!="D")
aerialB_noD = subset(aerialB, type!="D")
aerial_noD = subset(aerial, type!="D")

anova(lm(log10(spores+1)~ factor(spacing) + factor(dist.from.dis) + type,  
         data=aerialA_noD))

anova(lm(log10(spores+1)~ factor(spacing) + factor(dist.from.dis)*type,  
         data=aerialA_noD))

anova(lm(log10(spores+1)~ factor(spacing) + factor(dist.from.dis) + type,  
         data=aerialB_noD))

anova(lm(log10(spores+1)~ factor(spacing) + factor(dist.from.dis)*type,  
         data=aerialB_noD))

#table S9 - both experiments in one model, consistent results with analyzing separately
anova(lm(log10(spores+1)~ factor(spacing) + Experiment + factor(dist.from.dis)  + type,  
         data=aerial_noD))

anova(lm(log10(spores+1)~ factor(spacing) + Experiment + factor(dist.from.dis)*type,  
         data=aerial_noD))


##################
#### TABLE S10 ####

aerialA_Donly = subset(aerialA, type=="D")
aerialB_Donly = subset(aerialB, type=="D")
aerial_Donly = subset(aerial, type=="D")

anova(lm(log10(spores+1)~ spacing,  
         data=aerialA_Donly))

anova(lm(log10(spores+1)~ nsources,  
         data=aerialA_Donly))

anova(lm(log10(spores+1)~ spacing,  
         data=aerialB_Donly))

anova(lm(log10(spores+1)~ nsources,  
         data=aerialB_Donly))

# table s10 - both experiments combined, consistent with analyses separately
anova(lm(log10(spores+1)~ spacing + Experiment,  
         data=aerial_Donly))

anova(lm(log10(spores+1)~ nsources + Experiment,  
         data=aerial_Donly))


#### FIGURE 4 ####
#Mean aerial spore deposition at different X positions in the floral arrays.
ggplot(data=aerial, aes(x=xx, y=log10(spores+1))) + 
  stat_summary(aes(col=type), fun.data = "mean_se") + xlab("X position (m)") + ylab("Log10(Mean spore count+1)") +
  facet_wrap(~factor(spacing, levels=c(2,1,0.5,0.33,0.2), 
                     labels=c("Spacing=2m", "Spacing=1m","Spacing=0.5m","Spacing=0.33m","Spacing=0.2m")), nrow=1) +
  scale_color_manual(values=c("black","#78c679","#f768a1"), labels=c("Diseased","None","Healthy"), name="Surrounding Flowers")

#mean spores at each distance
aerial_sp_meandist = aerial %>% group_by(spacing, Experiment, dist.from.dis, type, nsources) %>%
  summarize(mean_spores = mean(spores), sterrlog = std.error(log10(spores+1)), n_samples = length(spores))


aer_mean_noD = aerial_sp_meandist[aerial_sp_meandist$type!="D",]
aer_mean_noD$spacing = as.numeric(aer_mean_noD$spacing)
aer_mean_noD = aer_mean_noD[aer_mean_noD$spacing<2,]
aer_mean_noD$log10mean_spores1 = log10(aer_mean_noD$mean_spores+1)

aer_mean_noD$spacing_f = factor(aer_mean_noD$spacing, levels=c(1,0.5,0.33,0.2))
aer_mean_noD$o.spacing = ordered(aer_mean_noD$spacing, levels=c(1,0.5,0.33,0.2))

#non-ordered factor for spacing for interpretation

#### TABLE 1b ####
aerialmeans_gam_log_wt = gam(log10mean_spores1~s(dist.from.dis, by= o.spacing, k=3) + dist.from.dis + Experiment + spacing_f, 
                             data=aer_mean_noD, method="REML", weights=n_samples)
gam.check(aerialmeans_gam_log_wt)
summary(aerialmeans_gam_log_wt)
anova(aerialmeans_gam_log_wt)


aerialmeans_gam_log_wt = gam(log10mean_spores1~s(dist.from.dis, by= o.spacing, k=3) + dist.from.dis + Experiment + o.spacing, 
                             data=aer_mean_noD, method="REML", weights=n_samples)
gam.check(aerialmeans_gam_log_wt)
summary(aerialmeans_gam_log_wt)
anova(aerialmeans_gam_log_wt)

plot_smooths(model = aerialmeans_gam_log_wt, series = dist.from.dis, comparison = Experiment, facet_terms = o.spacing) + 
  geom_errorbar(data= aer_mean_noD, width=0, alpha=1,
                aes(col=Experiment, x=dist.from.dis, ymin= log10mean_spores1-sterrlog, ymax= log10mean_spores1+sterrlog)) +
  geom_point(data= aer_mean_noD, alpha=1, size=2.5, col="black",
             aes(x=dist.from.dis, y= log10(mean_spores+1), 
                 shape=Experiment, fill=Experiment)) +
  facet_wrap(~ factor(o.spacing,labels=c("Spacing=1m","Spacing=0.5m","Spacing=0.33m","Spacing=0.2m")),
             nrow=1) +
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("A", "B")) +
  scale_fill_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_color_manual(values=c("#41b6c4","#fc8d59"), name="Replicate", labels=c("A", "B")) +
  scale_linetype_manual(values=c("solid","dashed"), name="Replicate", labels=c("A", "B")) + 
  theme_bw() + ylab("Log10(Mean spores + 1)") + xlab("Distance from nearest diseased plant (m)")




#### Aerial vs floral spore deposition against number of source plants ####

fl_subset_mean2 = fl_subset_mean[,c(1:5,7)]

#mean of spores on either side (grouping "type" E and H)

aer_subset_mean2 = aerial_sp_meandist %>% filter(type!="D") %>%
  group_by(spacing, Experiment, dist.from.dis, nsources) %>%
  summarize(mean_spores = mean(mean_spores), n_samples=sum(n_samples))

aer_subset_mean2$Rep=aer_subset_mean2$Experiment
fl_subset_mean2$trans_mode = "Floral"
aer_subset_mean2$trans_mode = "Aerial"

both_types_subset = rbind(fl_subset_mean2,aer_subset_mean2)
tail(both_types_subset)
 
both_types_subset$Rep = factor(both_types_subset$Rep)
both_types_subset$trans_mode = factor(both_types_subset$trans_mode)
both_types_subset$o.trans_mode = ordered(both_types_subset$trans_mode)

both_types_subset_2m = both_types_subset %>% filter(dist.from.dis==2 | dist.from.dis==2.0004 | dist.from.dis==2.04 | dist.from.dis==1.8)

both_types_subset_2m$o.Rep = ordered(both_types_subset_2m$Rep)

#### TABLE S11 now ####
nsrc_gam =gam(log10(mean_spores+1) ~ s(nsources, by=o.trans_mode, k=5) + s(nsources, by=Rep, k=5) + trans_mode + Rep, data=both_types_subset_2m, weights = n_samples, method="REML")
summary(nsrc_gam)
anova(nsrc_gam)

#### FIGURE 5 now ####
#The rate of spore deposition on each floral or aerial target as it varies 
#with the number of source plants at different densities. Each point represents 
#the mean spore deposition on targets at 2m distance from the source of disease.

#remove ordering of factors for plotting
nsrc_gam =gam(log10(mean_spores+1) ~ s(nsources, by=trans_mode, k=5) + s(nsources, by=Rep, k=5) + trans_mode + Rep, data=both_types_subset_2m, weights = n_samples, method="REML")

plot_smooths(model = nsrc_gam, series = nsources, comparison = trans_mode, facet_terms = Rep) + 
  geom_point(data= both_types_subset_2m, col="black", alpha=0.7, size=2.5, pch=21,
             aes(x=nsources, y= log10(mean_spores+1), fill=trans_mode)) +
  facet_wrap(~factor(Rep, labels=c("Replicate A","Replicate B")), scales="free") +
  xlab("Number of source diseased plants") + 
  ylab("Log10(Mean spore count + 1)") + 
  scale_color_manual(values=c("#4d9221","#c51b7d"), name="Deposition onto", labels=c("Spore-Traps", "Flowers")) +
  scale_fill_manual(values=c("#4d9221","#c51b7d"), name="Deposition onto", labels=c("Spore-Traps", "Flowers"))  +
  scale_linetype_manual(values=c("solid","dashed"), name="Deposition onto", labels=c("Spore-Traps", "Flowers")) +
  theme_bw()


#### Control plots ####

cont_spores = read.csv("control_spore_data.csv", header=T)

#now get the subset of points and take means and SE at 0, 2, 4, 6 m 
cont_spores_meanpos = cont_spores %>% filter(dist.from.dis==0 | dist.from.dis==2 | dist.from.dis==4 | dist.from.dis==6) %>%
  group_by(spacing, dist.from.dis) %>% summarize(mean_spores = mean(spores_mm2_new), sterrlog = std.error(log10(spores_mm2_new+1)), n_samples = length(spores_mm2_new))

subs_cont_spores_meanpos=cont_spores_meanpos
#subs_cont_spores_meanpos = subset(cont_spores_meanpos, spacing<4)
subs_cont_spores_meanpos$spacing_f = factor(subs_cont_spores_meanpos$spacing, levels=c(4,1,0.5,0.33))
subs_cont_spores_meanpos$o.spacing = ordered(subs_cont_spores_meanpos$spacing_f, levels=c(4,1,0.5,0.33))

hist(cont_spores$spores_mm2_new, xlim=c(0,200), breaks=50)
max(cont_spores$spores_mm2_new)
cont_spores[cont_spores$spores_mm2_new>50,]
mean(log10(cont_spores$spores_mm2_new[cont_spores$spacing==4]+1))
mean(cont_spores$spores_mm2_new[cont_spores$spacing==4])
std.error(cont_spores$spores_mm2_new[cont_spores$spacing==4])

mean(cont_spores$spores_mm2_new[cont_spores$spacing==0.33])
std.error(cont_spores$spores_mm2_new[cont_spores$spacing==0.33])

cont_spores$spacing_f = factor(cont_spores$spacing, levels=c(4,1,0.5,0.33))
cont_spores$o.spacing = ordered(cont_spores$spacing_f, levels=c(4,1,0.5,0.33))

#just the means at a subset of locations (0, 2, 4, 6 m; or just the middle location for spacing=4m)
#### TABLE S12 ####
#lm of mean spores deposited at 
#specific positions in the absence of diseased source plants 

cont_lm = lm(log10(mean_spores+1)~ spacing + dist.from.dis, 
   data=subs_cont_spores_meanpos)
summary(cont_lm)
Anova(cont_lm)
  
#### Figure S6 ####
#Floral spore deposition in the absence of diseased source plants. 
cont_spores_meanposall = cont_spores %>% filter(dist.from.dis==0 | dist.from.dis==2 | dist.from.dis==4 | dist.from.dis==6 | spacing==4) %>%
  group_by(spacing, dist.from.dis) %>% summarize(mean_spores = mean(spores_mm2_new), sterrlog = std.error(log10(spores_mm2_new+1)), n_samples = length(spores_mm2_new))

ggplot(aes(y=log10(mean_spores+1), x=dist.from.dis), data=cont_spores_meanposall) + 
  facet_wrap(~factor(spacing, levels=c(4, 1, 0.5, 0.33),
                     labels=c("Spacing=4m+","Spacing=1m","Spacing=0.5m","Spacing=0.33m")), nrow=1) +
  geom_errorbar(col="gray20", aes(x=dist.from.dis, ymin= log10(mean_spores+1)-sterrlog, ymax= log10(mean_spores+1)+sterrlog), width=0) + 
  geom_point(col="gray20", fill="gray", pch=22, size=2.5) + 
  theme_bw() + xlab("Distance (m)") + ylab("Log10(Mean spore count+1)") +
  theme(legend.position = "none")


#### Fan design data ####

# data from the fan design experiment
fandata = read.csv("fan_data.csv", header=T)

fandata_means = fandata %>% group_by(min_dist_dis, Experiment_replicate) %>% summarize(mean_count= mean(spores_mm2_new), n_samples=length(spores_mm2_new))

fan_gam1 = lm(log10(spores_mm2_new+1)~ s(min_dist_dis, k=9), 
                data=fandata, method="REML")
summary(fan_gam1)
gam.check(fan_gam1)
anova(fan_gam1)

#### Figure S7 #### 
#Floral spore deposition at increasing distances (m) from the nearest diseased flower in the fan designs

ggplot(data=fandata_means, aes(x=min_dist_dis, y=mean_count)) +
  #geom_smooth(aes(y=prediction), method="lm", formula=y~x, col="gray40") +
  geom_point( aes(pch=Experiment_replicate, fill=Experiment_replicate), size=2.5, alpha=0.7) + theme_bw() + xlab("Minimum distance from disease (m)") +
  ylab("Mean spore count") +
  scale_shape_manual(values=c(21,23),name="Replicate", labels=c("a", "b")) +
  scale_fill_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  scale_color_manual(values=c("#78c679","#9e9ac8"), name="Replicate", labels=c("a", "b")) +
  theme_bw() + ylab("Mean spore count") + xlab("Distance from disease (m)")+ xlim(0,10)


max(fandata_means$min_dist_dis)

