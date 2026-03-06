# This script runs group-wise comparison between a subset of VABS-3 matched TD and ASD (high adaptive)
# produce figure corresponding to supplementary figure 9

rm(list=ls())

# Load the required libraries
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(effectsize)
library(patchwork)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(Hmisc)

### Load data and parse table
# load the data
load("../preprocessed_data.RData")
load("../params_expdata_svnarrow.RData")

# reformat GEOxxxx_y/GEOxxxx as xxxx
alldata$gid <- as.integer(sub("^GEO(\\d+).*", "\\1", alldata$subID))

# create a binary group variable
alldata[, binary_group := fifelse(grepl("ASD", group, ignore.case = TRUE), "ASD", "TD")]
cat("allData N:",length(unique(alldata$subID)))

# create a verbal group variable
alldata[, group_verbal := fcase(
  is.na(vinabcstd), "flag",
  binary_group == "ASD" & vinabcstd >= 75, "ASD (high adaptive)",
  binary_group == "TD" & vinabcstd >= 75, "TD",
  binary_group == "ASD" & vinabcstd < 75, "ASD (low adaptive)",
  binary_group == "TD" & vinabcstd < 75, "flag_low_td",
  default = "Uncategorized"
)]

# create a srs normal group variable
alldata[, srs_normal := fifelse(
  (binary_group == "ASD" & srs2total > 20) | (binary_group == "TD" & srs2total < 120),1,0)]
# alldata[, table(srs_normal)]

# create a above criteria group variable
alldata[, above_crit := fifelse(
  (nTrials == 200 & perf > .616),1,0)]

# data summary
summary <- unique(alldata[, .(username, subID, game_version, source, binary_group, group_verbal, srs_normal, gender, race, age, ethnicity, vinabcstd, srs2total, bistotal,aasp_low_reg_raw,
                              aasp_sen_seek_raw,aasp_sen_ses_raw,aasp_sen_avoid_raw, perf, trainperf, nTrials, above_crit)])
cat("summary N:",length(unique(summary$subID)))

df = summary[game_version == "ft" & srs_normal==1 & group_verbal != "flag" & group_verbal != "flag_low_td" & above_crit==1]
df$group_verbal <- factor(as.character(df$group_verbal), levels = c("TD","ASD (high adaptive)","ASD (low adaptive)"),ordered = TRUE)

### Filter df for TD and ASD matched on VABS-3 scores:
df_overlap<-df%>%filter(vinabcstd>85 & vinabcstd<105)
cat("Supplementary figure df N:",length(unique(df$subID)))

df$is_subsample <- df$subID%in%unique(df_overlap$subID)

# supplementary figure p1: show range of ptpts
p1 <- ggplot(df,aes(x=vinabcstd,y=perf*100))+
  geom_point(aes(color=group_verbal,alpha = is_subsample+0.1))+
  ylab("%Correct")+
  xlab("VABS-3")+
  scale_color_manual(values = c("blue","orange","red"))+
  geom_vline(xintercept = 85,linetype="dashed",size=.5)+
  geom_vline(xintercept = 105,linetype="dashed",size=.5)+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


### Estimate psychometric slopes
# Fit a GLM model to the choice data for each subject
df_overlap_glm <- alldata[subID%in%unique(df_overlap$subID) & trial>10 & srs_normal == 1 & above_crit == 1 & game_version == "ft",c("username","subID","choice","dflash","winstay","loseswitch","binary_group","group_verbal")]
glm_id <- unique(df_overlap_glm$subID)

# Function to fit a GLM model to the choice data for a single subject
fit_subject_glm <- function(subject_data) {
  glm(choice ~ dflash + winstay + loseswitch, 
      data = subject_data, 
      family = binomial(link = "logit"))
}
# fit the model for each subject
model_fits <- df_overlap_glm %>%
  group_by(subID) %>%
  nest() %>%
  mutate(model = map(data, fit_subject_glm))

# get the parameter estimates for each subject
param_estimates_glm <- model_fits %>%
  mutate(params = map(model, broom::tidy)) %>%
  unnest(params) %>%
  select(subID, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)

# rename columnames
colnames(param_estimates_glm) <- c("subID","Sidebias","Slope","WS","LS")

# Merge with subject group information
param_estimates_glm <- param_estimates_glm %>% right_join(df_overlap, by = "subID")


### Extract fitted HBM parameters:
asdhbmfit <- rstan::extract(fitasd)
tdhbmfit <- rstan::extract(fittd)

params_df <- data.frame(group=rep(c("ASD","TD"),each=length(asdhbmfit$mu0)),
                        mu0=c(asdhbmfit$mu0,tdhbmfit$mu0),
                        mu1=c(asdhbmfit$mu1,tdhbmfit$mu1),
                        sig0=c(asdhbmfit$sig0,tdhbmfit$sig0),
                        sig1=c(asdhbmfit$sig1,tdhbmfit$sig1),
                        k0=c(asdhbmfit$k0,tdhbmfit$k0))
df = alldata[srs_normal == 1 & above_crit == 1 & game_version == "ft"]

df1 = data.frame("k1"=colMeans(asdhbmfit$k1),
                 "b1"=colMeans(asdhbmfit$b1),
                 "b2"=colMeans(asdhbmfit$b2),
                 "bs"=colMeans(asdhbmfit$bs),
                 "gr"="ASD")

df2 = data.frame("k1"=colMeans(tdhbmfit$k1),
                 "b1"=colMeans(tdhbmfit$b1),
                 "b2"=colMeans(tdhbmfit$b2),
                 "bs"=colMeans(tdhbmfit$bs),
                 "gr"="TD")

dt1 = df[binary_group=="ASD"]
dt2 = df[binary_group=="TD"]
dt1$subj_id<-as.numeric(factor(dt1$username, 
                               levels=unique(dt1$username)))

dt2$subj_id<-as.numeric(factor(dt2$username, 
                               levels=unique(dt2$username)))

df1$subj_id = as.numeric(c(1:length(unique(dt1$subj_id))))
df2$subj_id = as.numeric(c(1:length(unique(dt2$subj_id))))

dt1 = merge(dt1,df1,by="subj_id")
dt2 = merge(dt2,df2,by="subj_id")

dt = rbind(dt1,dt2)
dt = unique(dt[,.(username,k1,b1,b2,bs,binary_group)])

rm(dt1,dt2,df1,df2,df)

param_estimates_HBM <- dt

# HBM was fitted to the first attempt of the fixed-duration task
summary_ID = summary[game_version=="ft",c("username","binary_group","subID")]
# merge HBM parameters with demographic and survey info, matched by mongodb username
param_estimates_HBM = merge(summary_ID,param_estimates_HBM,by=c("username","binary_group"),all=TRUE)

# Merge with subject group information
param_glm_HBM <- param_estimates_glm %>% left_join(param_estimates_HBM, by = c("subID","binary_group","username"))

param_estimates <- as.data.table(param_glm_HBM)

param_estimates[, group_matched := fcase(
  group_verbal =="ASD (high adaptive)", "ASD (matched)",
  default = "TD"
)]

param_estimates$group_matched <- factor(as.character(param_estimates$group_matched), 
                                        levels = c("TD","ASD (matched)"),ordered = TRUE)


# VABS overlapped cohort analysis

# VABS-3 scores:
p2 <- ggplot(param_estimates,aes(x=group_matched,y=vinabcstd))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("VABS-3")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
        stat_compare_means(method = "t.test", paired = FALSE, 
                           method.args=c(var.equal=FALSE),
                           comparisons = list(c("TD", "ASD (matched)")), 
                           label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 
# run test - VABS-3
TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$vinabcstd
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$vinabcstd
sprintf("ASD N=%s",length(ASD_overlap))

# wilcoxon
wilcox.test(x=TD_overlap,
            y=ASD_overlap,
            paired = FALSE,exact=FALSE)
t.test(TD_overlap, ASD_overlap, var.equal = TRUE)

# Overall performance
p3 <- ggplot(param_estimates,aes(x=group_matched,y=perf*100))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("%Correct")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  geom_hline(yintercept = 62,linetype="dashed",size=.5)+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(50,100)+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) 

# run test - overall performance
TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$perf
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$perf
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)

# Slope
p4 <- ggplot(param_estimates,aes(x=group_matched,y=Slope))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("Psychometric Slope")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  geom_hline(yintercept = 0,linetype="dashed",size=.5)+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 

TD_overlap <- param_estimates[group_matched=="TD" & !is.na(Slope),Slope]
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(Slope),Slope]
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)


# HBM parameter k1 - Integration Noise
p5 <- ggplot(param_estimates,aes(x=group_matched,y=k1))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("Integration noise")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  geom_hline(yintercept = 0,linetype="dashed",size=.5)+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 
# test HBM k1 group difference between TD matched ASD and TD
ASD_overlap <- param_estimates[group_matched=="ASD (matched)",]$k1
TD_overlap <- param_estimates[group_matched=="TD",]$k1
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)

# Additional survey measures analyses:
# SRS-2
p6 <- ggplot(param_estimates,aes(x=group_matched,y=srs2total))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("SRS-2 total")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 
# test SRS-2
TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$srs2total
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$srs2total
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)

# BIS
p7 <- ggplot(param_estimates,aes(x=group_matched,y=bistotal))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("BIS total")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 
# run test - bis total score
TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$bistotal
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$bistotal
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)

# AASP
p8 <- ggplot(param_estimates,aes(x=group_matched,y=aasp_low_reg_raw))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("AASP (low registration)")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 

TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$aasp_low_reg_raw
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$aasp_low_reg_raw
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)


p9 <- ggplot(param_estimates,aes(x=group_matched,y=aasp_sen_avoid_raw))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("AASP (sensation avoiding)")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 

TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$aasp_sen_avoid_raw
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$aasp_sen_avoid_raw
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)

p10 <- ggplot(param_estimates,aes(x=group_matched,y=aasp_sen_seek_raw))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("AASP (sensation seeking)")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 

TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$aasp_sen_seek_raw
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$aasp_sen_seek_raw
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)

p11 <- ggplot(param_estimates,aes(x=group_matched,y=aasp_sen_ses_raw))+
  geom_beeswarm(aes(color=group_matched))+
  ylab("AASP (sensory sensitivity)")+
  xlab("Group")+
  scale_color_manual(values = c("blue","orange","red"))+
  theme_bw()+
  labs(color="")+
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(method = "t.test", paired = FALSE, 
                     method.args=c(var.equal=FALSE),
                     comparisons = list(c("TD", "ASD (matched)")), 
                     label = "p.signif")+
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) 

TD_overlap <- param_estimates[group_matched=="TD" & !is.na(perf),]$aasp_sen_ses_raw
sprintf("TD N=%s",length(TD_overlap))
ASD_overlap <- param_estimates[group_matched=="ASD (matched)" & !is.na(perf),]$aasp_sen_ses_raw
sprintf("ASD N=%s",length(ASD_overlap))
# wilcoxon test
wilcox.test(x=TD_overlap, y=ASD_overlap, paired = FALSE, exact = FALSE)
# Welch two sample t-test
t.test(TD_overlap, ASD_overlap, var.equal = FALSE)
# effect size:
effectsize::cohens_d(TD_overlap,ASD_overlap)

fig <- (p1 + p2)/(p3 + p4 + p5)/(p6 | p7)/(p8 | p9 | p10 | p11) + plot_annotation(tag_levels = 'A')
ggsave(filename = "./supplementary_matchedVABS-3_Welch_t-test.pdf", plot = fig, width = 10, height = 12,dpi = 'retina')
