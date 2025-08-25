############################################################
## Project: Statistical Modeling of Decision-Making & fMRI — Original Pipeline (Organized)
## File: decision_fMRI_pipeline.R
## Author: Anthony Romyn
##
## Dependencies:
##   lme4, lmerTest, devtools, nifti.io, tkmisc, nifti.cluster, mRi,
##   ggplot2, moments, ggeffects, lavaan
############################################################

## ---------------------------
## 0) Libraries
## ---------------------------
library(lme4)
library(lmerTest)
library(devtools)
# (Run these install_github lines only if needed)
# devtools::install_github("TKoscik/mRi", auth_token="<redacted>")
# devtools::install_github("TKoscik/mRi", auth_token="<redacted>")
devtools::install_github("TKoscik/nifti.io")
devtools::install_github("TKoscik/tkmisc")
devtools::install_github("TKoscik/nifti.cluster")

library(nifti.io)
library(tkmisc)
library(nifti.cluster)
library(mRi)
library(ggplot2)
library(moments)

## ---------------------------
## 1) Load data
## ---------------------------
df <- read.csv("./dualspin/stbetas/dualSpin.csv", header = TRUE)
ptplist <- read.csv("./dualspin/stbetas/ptplist.txt", header = FALSE)
colnames(ptplist) <- "ptp"

NAcc_cluster <- read.csv(
  "./dualspin/stbetas/output/decision_3_anteriorMFG_lateralOFC/cluster/vol18clusters_p001cl6.csv",
  header = TRUE
)

mainEffect_clusters <- read.csv(
  "./dualspin/stbetas/output/decision_1_2/cluster/vol18clusters_p001sz117cl6.csv",
  header = TRUE
)

vmpfc_selfAm <- read.csv(
  "./dualspin/stbetas/neuralEmpathyMeasures/clusteredAtp001/allPtps_clust_0003_selfAm.1D",
  header = FALSE
)
vmpfc_otherAm <- read.csv(
  "./dualspin/stbetas/neuralEmpathyMeasures/clusteredAtp001/allPtps_clust_0003_otherAm.1D",
  header = FALSE
)
colnames(vmpfc_selfAm) <- "vmpfc_selfAm"
colnames(vmpfc_otherAm) <- "vmpfc_otherAm"

## ---------------------------
## 2) Merge & centering
## ---------------------------
df2 <- cbind(ptplist, vmpfc_selfAm, vmpfc_otherAm)
df2$vmpfc_selfAm_c  <- scale(df2$vmpfc_selfAm,  scale = TRUE, center = TRUE)
df2$vmpfc_otherAm_c <- scale(df2$vmpfc_otherAm, scale = TRUE, center = TRUE)

df3 <- merge(df, df2, by = "ptp")
df3 <- cbind(df3, NAcc_cluster, mainEffect_clusters)
df  <- df3
df  <- df[which(is.na(df$takepass_contrast) != 1), ]
colnames(df)

## ---------------------------
## 3) Decision model: GLMM (binomial)
## ---------------------------

# Remove outliers and decompose variable
ex1 <- outliers.by.subject(df, response = "anteriorMFG_lateralOFC_r", subject = "ptp", trim = 3)
df  <- decompose.var(ex1$data, var.name = "anteriorMFG_lateralOFC_r", grouping = c("grand", "ptp", "Run"))

tf2 <- df
tf2$vmpfc_selfAm_centered  <- scale(tf2$vmpfc_selfAm,  scale = TRUE, center = TRUE)
tf2$vmpfc_otherAm_centered <- scale(tf2$vmpfc_otherAm, scale = TRUE, center = TRUE)

mA2 <- glmer(
  Decision ~ (selfEV + otherEV) * (anteriorMFG_lateralOFC_r.trial.c) * (vmpfc_selfAm_centered + vmpfc_otherAm_centered) +
    (1 + selfEV + otherEV + anteriorMFG_lateralOFC_r.trial.c | ptp/Run),
  data = tf2,
  family = binomial,
  control = glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap", optCtrl = list(maxfun = 1000000))
)
summary(mA2)

# Remove outliers by standardized residuals and re-decompose
ex2 <- outliers.by.std.resid(model = mA2, data = tf2, trim = 3)
tf3 <- ex2$data
tf3$anteriorMFG_lateralOFC_r.ptp.c   <- NULL
tf3$anteriorMFG_lateralOFC_r.Run.c   <- NULL
tf3$anteriorMFG_lateralOFC_r.trial.c <- NULL
tf3 <- decompose.var(tf3, var.name = "anteriorMFG_lateralOFC_r", grouping = c("grand", "ptp", "Run"))

tf3$vmpfc_selfAm_centered  <- scale(tf3$vmpfc_selfAm,  scale = TRUE, center = TRUE)
tf3$vmpfc_otherAm_centered <- scale(tf3$vmpfc_otherAm, scale = TRUE, center = TRUE)

mA3 <- glmer(
  Decision ~ (selfEV + otherEV) * (anteriorMFG_lateralOFC_r.trial.c) * (vmpfc_selfAm_centered + vmpfc_otherAm_centered) +
    (1 + selfEV + otherEV + anteriorMFG_lateralOFC_r.trial.c | ptp/Run),
  data = tf3,
  family = binomial,
  control = glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap", optCtrl = list(maxfun = 1000000))
)
summary(mA3)

mdl_predDecision <- mA3

## ---------------------------
## 4) Decision model: Effects plot
## ---------------------------
library(ggeffects)
library(ggplot2)

sd(tf3$anteriorMFG_lateralOFC_r.trial.c)

mdl_predDecisionfit <- ggpredict(
  mdl_predDecision,
  terms = c("otherEV [all]", "anteriorMFG_lateralOFC_r.trial.c [-0.385,0,0.385]", "vmpfc_otherAm_centered [-1,1]")
)
levels(mdl_predDecisionfit$facet) <- c("-1 SD Empathy", "+1 SD Empathy")

ggplot(mdl_predDecisionfit, aes(x, y = predicted, group = group, colour = group, fill = group)) +
  facet_grid(~ facet) +
  geom_line() +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low), alpha = 0.3, linetype = 0, show.legend = FALSE) +
  scale_x_continuous(name = "other gamble expected value", breaks = seq(-5, 5, 5)) +
  theme(panel.grid.major = element_line(colour = "white")) +
  scale_y_continuous(name = "Decision", breaks = c(0.1, 0.5, 0.95), labels = paste0(c("Pass", "50%", "Take"))) +
  scale_color_manual(name = "VLPFC", labels = c("-1 SD", "0", "+1 SD"), values = c("#D55E00", "#0072B2", "#009E73")) +
  scale_fill_manual(labels = c("-1 SD", "0", "+1 SD"), values = c("#D55E00", "#0072B2", "#009E73")) +
  ggtitle("")

## ---------------------------
## 5) Predicting NAcc (LMM)
## ---------------------------

ex1 <- outliers.by.subject(df, response = "otherEV_dlpfc_neuralEmpathy_predNacc", subject = "ptp", trim = 3)
ex1 <- outliers.by.subject(ex1$data, response = "anteriorMFG_lateralOFC_r", subject = "ptp", trim = 3)
ex1$data <- decompose.var(ex1$data, var.name = "anteriorMFG_lateralOFC_r", grouping = c("grand", "ptp", "Run"))
tf2 <- ex1$data

mA2 <- lmer(
  otherEV_dlpfc_neuralEmpathy_predNacc ~ (selfEV + otherEV) * anteriorMFG_lateralOFC_r.trial.c * (vmpfc_selfAm_c + vmpfc_otherAm_c) +
    (1 + selfEV + otherEV + anteriorMFG_lateralOFC_r.trial.c | ptp/Run),
  data = tf2,
  control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap", optCtrl = list(maxfun = 1000000))
)
summary(mA2)

ex2 <- outliers.by.std.resid(model = mA2, data = tf2, trim = 3)
tf3 <- ex2$data
colnames(tf3)
tf3$anteriorMFG_lateralOFC_r.trial.c <- NULL
tf3$anteriorMFG_lateralOFC_r.ptp.c   <- NULL
tf3$anteriorMFG_lateralOFC_r.Run.c   <- NULL
colnames(tf3)
tf3 <- decompose.var(tf3, var.name = "anteriorMFG_lateralOFC_r", grouping = c("grand", "ptp", "Run"))
colnames(tf3)
tf3$vmpfc_selfAm_c  <- scale(tf3$vmpfc_selfAm,  scale = TRUE, center = TRUE)
tf3$vmpfc_otherAm_c <- scale(tf3$vmpfc_otherAm, scale = TRUE, center = TRUE)

mA3 <- lmer(
  otherEV_dlpfc_neuralEmpathy_predNacc ~ (selfEV + otherEV) * anteriorMFG_lateralOFC_r.trial.c * (vmpfc_selfAm_c + vmpfc_otherAm_c) +
    (1 + selfEV + otherEV + anteriorMFG_lateralOFC_r.trial.c | ptp/Run),
  data = tf3,
  control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap", optCtrl = list(maxfun = 1000000))
)
summary(mA3)

sd(tf3$anteriorMFG_lateralOFC_r.trial.c)
sd(tf3$vmpfc_otherAm_c)

fit <- ggpredict(
  mA3,
  terms = c("otherEV [all]", "anteriorMFG_lateralOFC_r.trial.c [-0.385,0,0.385]", "vmpfc_otherAm_c [-1,1]")
)
levels(fit$facet) <- c("-1 SD Empathy", "+1 SD Empathy")

ggplot(fit, aes(x, y = predicted, group = group, colour = group, fill = group)) +
  facet_grid(~ facet) +
  geom_line() +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low), alpha = 0.3, linetype = 0, show.legend = FALSE) +
  scale_x_continuous(name = "other gamble expected value", breaks = seq(-5, 5, 5)) +
  theme(panel.grid.major = element_line(colour = "white")) +
  scale_y_continuous(name = "NAcc") +
  scale_color_manual(name = "VLPFC", labels = c("-1 SD", "0", "+1 SD"), values = c("#D55E00", "#0072B2", "#009E73")) +
  scale_fill_manual(labels = c("-1 SD", "0", "+1 SD"), values = c("#D55E00", "#0072B2", "#009E73")) +
  ggtitle("")

## ---------------------------
## 6) Double-mean centering approach for SEM (lavaan)
## ---------------------------
library(lavaan)

# Make sure all measured vars are mean-centered
tf2$NAcc_cluster_c <- scale(tf2$otherEV_dlpfc_neuralEmpathy_predNacc, scale = TRUE, center = TRUE)[, ]
tf2$vmpfc_selfAm_c_centered   <- scale(tf2$vmpfc_selfAm_c,   scale = FALSE, center = TRUE)
tf2$vmpfc_otherAm_c_centered  <- scale(tf2$vmpfc_otherAm_c,  scale = FALSE, center = TRUE)
tf2$selfEV_centered           <- scale(tf2$selfEV,           scale = FALSE, center = TRUE)
tf2$otherEV_centered          <- scale(tf2$otherEV,          scale = FALSE, center = TRUE)
tf2$anteriorMFG_lateralOFC_r.trial.c_centered <- scale(tf2$anteriorMFG_lateralOFC_r.trial.c, scale = FALSE, center = TRUE)

# Interactions
tf2$selfEVByDLPFC         <- tf2$selfEV_centered  * tf2$anteriorMFG_lateralOFC_r.trial.c_centered
tf2$otherEVByDLPFC        <- tf2$otherEV_centered * tf2$anteriorMFG_lateralOFC_r.trial.c_centered
tf2$selfEVByvmpfc_self    <- tf2$selfEV_centered  * tf2$vmpfc_selfAm_c_centered
tf2$selfEVByvmpfc_other   <- tf2$selfEV_centered  * tf2$vmpfc_otherAm_c_centered
tf2$otherEVByvmpfc_self   <- tf2$otherEV_centered * tf2$vmpfc_selfAm_c_centered
tf2$otherEVByvmpfc_other  <- tf2$otherEV_centered * tf2$vmpfc_otherAm_c_centered

tf2$selfEVByDLPFCByvmpfc_self  <- tf2$selfEV_centered  * tf2$anteriorMFG_lateralOFC_r.trial.c_centered * tf2$vmpfc_selfAm_c_centered
tf2$selfEVByDLPFCByvmpfc_other <- tf2$selfEV_centered  * tf2$anteriorMFG_lateralOFC_r.trial.c_centered * tf2$vmpfc_otherAm_c_centered
tf2$otherEVByDLPFCByvmpfc_self <- tf2$otherEV_centered * tf2$anteriorMFG_lateralOFC_r.trial.c_centered * tf2$vmpfc_selfAm_c_centered
tf2$otherEVByDLPFCByvmpfc_other<- tf2$otherEV_centered * tf2$anteriorMFG_lateralOFC_r.trial.c_centered * tf2$vmpfc_otherAm_c_centered

# Mean-center each interaction
tf2$selfEVByDLPFC         <- scale(tf2$selfEVByDLPFC,         scale = FALSE, center = TRUE)
tf2$otherEVByDLPFC        <- scale(tf2$otherEVByDLPFC,        scale = FALSE, center = TRUE)
tf2$selfEVByvmpfc_self    <- scale(tf2$selfEVByvmpfc_self,    scale = FALSE, center = TRUE)
tf2$selfEVByvmpfc_other   <- scale(tf2$selfEVByvmpfc_other,   scale = FALSE, center = TRUE)
tf2$otherEVByvmpfc_self   <- scale(tf2$otherEVByvmpfc_self,   scale = FALSE, center = TRUE)
tf2$otherEVByvmpfc_other  <- scale(tf2$otherEVByvmpfc_other,  scale = FALSE, center = TRUE)

tf2$selfEVByDLPFCByvmpfc_self   <- scale(tf2$selfEVByDLPFCByvmpfc_self,   scale = FALSE, center = TRUE)
tf2$selfEVByDLPFCByvmpfc_other  <- scale(tf2$selfEVByDLPFCByvmpfc_other,  scale = FALSE, center = TRUE)
tf2$otherEVByDLPFCByvmpfc_self  <- scale(tf2$otherEVByDLPFCByvmpfc_self,  scale = FALSE, center = TRUE)
tf2$otherEVByDLPFCByvmpfc_other <- scale(tf2$otherEVByDLPFCByvmpfc_other, scale = FALSE, center = TRUE)

# SEM model
mod.med.lavaan2 <- '
Decision ~ cdash1*selfEV_centered + cdash2*otherEV_centered + cdash3*anteriorMFG_lateralOFC_r.trial.c_centered + cdash4*vmpfc_selfAm_c_centered + cdash5*vmpfc_otherAm_c_centered + cdash6*selfEVByDLPFC + cdash7*selfEVByvmpfc_self + cdash8*selfEVByvmpfc_other + cdash9*otherEVByDLPFC + cdash10*otherEVByvmpfc_self + cdash11*otherEVByvmpfc_other + cdash12*selfEVByDLPFCByvmpfc_self + cdash13*selfEVByDLPFCByvmpfc_other + cdash14*otherEVByDLPFCByvmpfc_self + cdash15*otherEVByDLPFCByvmpfc_other

NAcc_cluster_c ~ a1*selfEV_centered + a2*otherEV_centered + a3*anteriorMFG_lateralOFC_r.trial.c_centered + a4*vmpfc_selfAm_c_centered + a5*vmpfc_otherAm_c_centered + a6*selfEVByDLPFC + a7*selfEVByvmpfc_self + a8*selfEVByvmpfc_other + a9*otherEVByDLPFC + a10*otherEVByvmpfc_self + a11*otherEVByvmpfc_other + a12*selfEVByDLPFCByvmpfc_self + a13*selfEVByDLPFCByvmpfc_other + a14*otherEVByDLPFCByvmpfc_self  + a15*otherEVByDLPFCByvmpfc_other

Decision ~ b1*NAcc_cluster_c

vmpfc_otherAm_c_centered ~ vmpfc_otherAm_c_centered.mean*1
vmpfc_otherAm_c_centered ~~ vmpfc_otherAm_c_centered.var*vmpfc_otherAm_c_centered
vmpfcSD:=sqrt(vmpfc_otherAm_c_centered.var)

indirect.SDbelow := (a9 + a15*(vmpfc_otherAm_c_centered.mean-(1*vmpfcSD)))*b1
indirect.SDabove := (a9 + a15*(vmpfc_otherAm_c_centered.mean+(1*vmpfcSD)))*b1

direct.SDbelow := cdash9 + cdash15*(vmpfc_otherAm_c_centered.mean-(1*vmpfcSD)) 
direct.SDabove := cdash9 + cdash15*(vmpfc_otherAm_c_centered.mean+(1*vmpfcSD))

total.SDbelow := direct.SDbelow + indirect.SDbelow
total.SDabove := direct.SDabove + indirect.SDabove

prop.mediated.SDbelow := indirect.SDbelow / total.SDbelow
prop.mediated.SDabove := indirect.SDabove / total.SDabove

index.mod.med := a15*b1
'

fit.mod.med.lavaan2 <- sem(mod.med.lavaan2, tf2, test = "bootstrap", bootstrap = 200)
summary(fit.mod.med.lavaan2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
