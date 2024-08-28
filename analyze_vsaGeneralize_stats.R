library(effectsize)
library(multcomp)
library(circular)
library(circglmbayes)

setwd("/Volumes/smng/experiments/vsaGeneralize/acousticdata")

## global (AVS) statistics
AVS <-read.csv("datatable_AVS.txt")
adapt <- AVS[ which(AVS$phase=='adapt'),] 
gen <- AVS[ which(AVS$phase=='gen'),] 

#report summary stats
c(mean(adapt$AVS),sd(adapt$AVS)/sqrt(length(adapt$AVS)),mean(gen$AVS),sd(gen$AVS)/sqrt(length(gen$AVS)))

#t-test for difference from 0 and effect sizes
t.test(adapt$AVS)
hedges_g(adapt$AVS)
t.test(gen$AVS)
hedges_g(gen$AVS)

#correlation between adaptation and generalization
cor.test(adapt$AVS,gen$AVS,method="spearman")

## vowel-specific statistics, only looking at changes from baseline
dat <-read.csv("datatable_vowels.txt")
train <- dat[which(dat$vowel=='iy' | dat$vowel=='ae' | dat$vowel=='uw' | dat$vowel=='aa'),]
gen <- dat[which(dat$vowel=='ih' | dat$vowel=='ey' | dat$vowel=='eh' | dat$vowel=='ah' | dat$vowel=='ow'),]
gen$vowel <- factor(gen$vowel)

#individual t-tests
iy <- dat[which(dat$vowel=='iy'),]
ih <- dat[which(dat$vowel=='ih'),]
ey <- dat[which(dat$vowel=='ey'),]
eh <- dat[which(dat$vowel=='eh'),]
ae <- dat[which(dat$vowel=='ae'),]
aa <- dat[which(dat$vowel=='aa'),]
ah <- dat[which(dat$vowel=='ah'),]
ow <- dat[which(dat$vowel=='ow'),]
uw <- dat[which(dat$vowel=='uw'),]

pval = c(1:9);
pval[1] <- t.test(iy$mag)$p.value
pval[2] = t.test(ih$mag)$p.value
pval[3] = t.test(ey$mag)$p.value
pval[4] = t.test(eh$mag)$p.value
pval[5] = t.test(ae$mag)$p.value
pval[6] = t.test(aa$mag)$p.value
pval[7] = t.test(ah$mag)$p.value
pval[8] = t.test(ow$mag)$p.value
pval[9] = t.test(uw$mag)$p.value

pval_adj = p.adjust(pval,method="fdr")

mean(dat$mag)
sd(dat$mag)/sqrt(length((dat$mag)))

#anova for magnitude of adaptation and effect of distance
mag.aov = aov(mag~vowel+Error(subj),dat)
summary(mag.aov)
eta_squared(mag.aov)
pairwise.t.test(dat$mag,dat$vowel,paired=TRUE,p.adjust.method = "fdr")
#lmer generates singular fits with random intercepts for participants
mag.dist = lm(mag~trainDistMag,gen)
anova(mag.dist)

#as above, for angle are done in Matlab as the assumptions for means testing in 'circular' package aren't met. 
#see calc_vsaGeneralize_circStats.mat

## vowel-specific statistics, including baseline and end phases
dat_phase <-read.csv("datatable_vowels_byphase.txt")

#anova for magnitude of adaptation and effect of distance
mag.aov = aov(mag~vowel*phase+Error(subj),dat_phase)
summary(mag.aov)
eta_squared(mag.aov)



#individual t-tests
iy <- dat_phase[which(dat_phase$vowel=='iy'),]
ih <- dat_phase[which(dat_phase$vowel=='ih'),]
ey <- dat_phase[which(dat_phase$vowel=='ey'),]
eh <- dat_phase[which(dat_phase$vowel=='eh'),]
ae <- dat_phase[which(dat_phase$vowel=='ae'),]
aa <- dat_phase[which(dat_phase$vowel=='aa'),]
ah <- dat_phase[which(dat_phase$vowel=='ah'),]
ow <- dat_phase[which(dat_phase$vowel=='ow'),]
uw <- dat_phase[which(dat_phase$vowel=='uw'),]

pval_phase = c(1:9);
pval_phase[1] = t.test(iy[which(iy$phase=='Base'),]$mag,iy[which(iy$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[2] = t.test(ih[which(ih$phase=='Base'),]$mag,ih[which(ih$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[3] = t.test(ey[which(ey$phase=='Base'),]$mag,ey[which(ey$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[4] = t.test(eh[which(eh$phase=='Base'),]$mag,eh[which(eh$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[5] = t.test(ae[which(ae$phase=='Base'),]$mag,ae[which(ae$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[6] = t.test(aa[which(aa$phase=='Base'),]$mag,aa[which(aa$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[7] = t.test(ah[which(ah$phase=='Base'),]$mag,ah[which(ah$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[8] = t.test(ow[which(ow$phase=='Base'),]$mag,ow[which(ow$phase=='End'),]$mag,paired = TRUE)$p.value
pval_phase[9] = t.test(uw[which(uw$phase=='Base'),]$mag,uw[which(uw$phase=='End'),]$mag,paired = TRUE)$p.value

pval_phase_adj = p.adjust(pval_phase,method="fdr")


## clear speech characteristics statistics, including baseline and end phases
dat_clear <-read.csv("datatable_clearSpeech.txt")

#anova for vowel duration
dur.aov = aov(dur~vowel*phase+Error(subj/(vowel*phase)),dat_clear)
summary(dur.aov)
eta_squared(mag.aov)

#anova for intensity
int.aov = aov(int~vowel*phase+Error(subj/(vowel*phase)),dat_clear)
summary(int.aov)
eta_squared(int.aov)

#anova for f0 Max
f0Max.aov = aov(f0Max~vowel*phase+Error(subj/(vowel*phase)),dat_clear)
summary(f0Max.aov)
eta_squared(f0Max.aov)

#anova for f0 Range
f0Range.aov = aov(f0Range~vowel*phase+Error(subj/(vowel*phase)),dat_clear)
summary(f0Range.aov)
eta_squared(f0Range.aov)






