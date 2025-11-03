setwd("...")

library(ggpubr)
library(car)
library(rstatix)
library(ARTool)
library(emmeans)

# Version 1: Study Forrest Brightness and NormDiff ####
lowlvl_visualfeature = read.csv('./_Visual_Features/LowVisualFeature.csv')
movie_features = read.csv("./01_Stimuli/AdditionalFiles/trial_info_VisualFeature.csv")

movie_features$hierarchy = factor(movie_features$hierarchy, levels = c("Shot","Scene"))
movie_features$duration = factor(movie_features$duration)


names(lowlvl_visualfeature)
movie_features2 =aggregate(cbind(lum, contrast, saturation, complexity) ~ ContLvl+Duration+ trial+place+time+location, data = lowlvl_visualfeature, mean)

movie_features2$ContLvl = factor(movie_features2$ContLvl, levels = c("Shot","Scene"))
movie_features2$Duration = factor(movie_features2$Duration, level = c("4s","12s","36s"))

# Luminance ####
ggboxplot(movie_features2, x = "Duration", y = "lum", fill = "ContLvl")
# Normalverteilung der residuen
nv <- lm(lum ~ ContLvl*Duration, data = movie_features2)
ggqqplot(residuals(nv))

# Varianzhomogenität
leveneTest(data = movie_features2, lum ~ ContLvl*Duration)
plot(nv,1)

# ARTools (non-parametric alternative)
art_model <- art(lum ~ ContLvl*Duration, data = movie_features2)
m.art.anova = anova(art_model)

m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova


emmeans(artlm(art_model, "ContLvl"), pairwise ~ ContLvl, adjust = "bonf")
emmeans(artlm(art_model, "Duration"), pairwise ~ Duration, adjust = "bonf")


# Contrast ####
ggboxplot(movie_features2, x = "Duration", y = "contrast", fill = "ContLvl")
# Normalverteilung der residuen
nv <- lm(contrast ~ ContLvl*Duration, data = movie_features2)
ggqqplot(residuals(nv)) # violated

# Varianzhomogenität
leveneTest(data = movie_features2, contrast ~ ContLvl*Duration)
plot(nv,1) # violated

# ARTools (non-parametric alternative)
art_model <- art(contrast ~ ContLvl*Duration, data = movie_features2)
m.art.anova = anova(art_model)

m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova

emmeans(artlm(art_model, "ContLvl"), pairwise ~ ContLvl, adjust = "bonf")
emmeans(artlm(art_model, "Duration"), pairwise ~ Duration, adjust = "bonf")


# saturation ####
ggboxplot(movie_features2, x = "Duration", y = "saturation", fill = "ContLvl")
# Normalverteilung der residuen
nv <- lm(saturation ~ ContLvl*Duration, data = movie_features2)
ggqqplot(residuals(nv)) # violated

# Varianzhomogenität
leveneTest(data = movie_features2, saturation ~ ContLvl*Duration)
plot(nv,1) # tendency to violated

art_model <- art(saturation ~ ContLvl*Duration, data = movie_features2)
m.art.anova = anova(art_model)

m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova

emmeans(artlm(art_model, "ContLvl"), pairwise ~ ContLvl, adjust = "bonf")
emmeans(artlm(art_model, "Duration"), pairwise ~ Duration, adjust = "bonf")


# complexity ####
ggboxplot(movie_features2, x = "Duration", y = "complexity", fill = "ContLvl")
# Normalverteilung der residuen
nv <- lm(complexity ~ ContLvl*Duration, data = movie_features2)
ggqqplot(residuals(nv)) 
# Varianzhomogenität
leveneTest(data = movie_features2, complexity ~ ContLvl*Duration)
plot(nv,1) # okayish

# ANOVA
art_model <- art(complexity ~ ContLvl*Duration, data = movie_features2)
m.art.anova = anova(art_model)

m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova

emmeans(artlm(art_model, "ContLvl"), pairwise ~ ContLvl, adjust = "bonf")
emmeans(artlm(art_model, "Duration"), pairwise ~ Duration, adjust = "bonf")

# motion energy ####
motion_energy = read.csv('./_Visual_Features/Motion_energy.csv')
motion_energy$ContLvl = factor(motion_energy$ContLvl, levels = c("Shot","Scene"))
motion_energy$Duration = factor(motion_energy$Duration, level = c("4s","12s","36s"))

ggboxplot(motion_energy, x = "Duration", y = "motion_energy", fill = "ContLvl")
# Normalverteilung der residuen
nv <- lm(motion_energy ~ ContLvl*Duration, data = motion_energy)
ggqqplot(residuals(nv))

# Varianzhomogenität
leveneTest(data = motion_energy, motion_energy ~ ContLvl*Duration)
plot(nv,1)

art_model <- art(motion_energy ~ ContLvl*Duration, data = motion_energy)
m.art.anova = anova(art_model)

m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova

emmeans(artlm(art_model, "ContLvl"), pairwise ~ ContLvl, adjust = "bonf")
emmeans(artlm(art_model, "Duration"), pairwise ~ Duration, adjust = "bonf")

# normdiff####
ggboxplot(movie_features, x = "duration", y = "normdiff", fill = "hierarchy")
# Normalverteilung der residuen
nv <- lm(normdiff ~ hierarchy*duration, data = movie_features)
ggqqplot(residuals(nv))

# Varianzhomogenität
leveneTest(data = movie_features, normdiff ~ hierarchy*duration)
plot(nv,1)

anova_test(data = movie_features, normdiff ~ hierarchy*duration, effect.size = "pes") # partial eta squared

art_model <- art(normdiff ~ hierarchy * duration, data = movie_features)
m.art.anova = anova(art_model)

m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova

emmeans(artlm(art_model, "hierarchy"), pairwise ~ hierarchy, adjust = "bonf")
emmeans(artlm(art_model, "duration"), pairwise ~ duration, adjust = "bonf")

## Combined ####
movie_features2 = movie_features2[order(movie_features2$trial),]
movie_features2$normdiff = movie_features$normdiff
movie_features2$motion_energy = motion_energy$motion_energy

library(dplyr)

# Define your dependent variables
dvs <- c("lum", "contrast", "saturation", "complexity", "normdiff", "motion_energy")

# Create an empty list to store results
results <- list()
pvals <- matrix(data = NA, nrow = length(dvs), ncol = 3)
DF = c()
DF_RES = c()
F_VALUES = c()
PART_ETA_SQ = c()
m = 0
for (dv in dvs) {
  m = m +1 
  formula <- as.formula(paste(dv, "~ ContLvl * Duration"))
  
  model <- art(formula, data = movie_features2)
  m.art.anova = anova(model)
  
  m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
  
  
  m.art.anova$DV <- dv  # tag variable name
  results[[dv]] <- m.art.anova
  
  # Extract p-values for correction later
  pvals[m, ] = m.art.anova$`Pr(>F)`
  DF = c(DF,m.art.anova$Df)
  DF_RES = c(DF_RES, m.art.anova$Df.res)
  F_VALUES = c(F_VALUES,m.art.anova$`F value`)
  PART_ETA_SQ = c(PART_ETA_SQ,m.art.anova$eta.sq.part)
  
  cat("\n==============================\n")
  cat("Dependent variable:", dv, "\n")
  print(art_aov)
}

# Optional: multiple-comparison correction
pvals_corrected = c()
pvals_corrected = cbind(pvals_corrected, p.adjust(pvals[,1], method = "fdr"))
pvals_corrected = cbind(pvals_corrected, p.adjust(pvals[,2], method = "fdr"))
pvals_corrected = cbind(pvals_corrected, p.adjust(pvals[,3], method = "fdr"))
pvals_corrected_shape = as.vector(t(pvals_corrected))


pvals_uncorrected = as.vector(t(pvals))


# Combine into a summary table
summary_df <- data.frame(
  DV = rep(dvs, each = 3),  # 3 terms per DV: ContLvl, Duration, ContLvl:Duration
  Term = rep(c("ContLvl", "Duration", "Interaction"), times = length(dvs)),
  f_value = round(F_VALUES,2),
  df1 = DF,
  df2 = DF_RES,
  p_value = round(pvals_uncorrected,3),
  p_value_FDR = round(pvals_corrected_shape,3),
  eta_sq_part = round(PART_ETA_SQ,2)
)

print(summary_df)

write.csv(summary_df, "./results_low_level_features.csv")


# Scene/Shot Position alogn movie runtime ####
ggplot(movie_features) +
  geom_segment(aes(x = frame_position, xend = frame_position, y = 0, yend = 1, 
                   colour = interaction(hierarchy, duration)), alpha = 1) +
  theme_minimal() + facet_wrap(~hierarchy*duration, scales = "free")

ggboxplot(movie_features, x = "duration", y = "frame_position", fill = "hierarchy")
# Normalverteilung der residuen
nv <- lm(frame_position ~ hierarchy*duration, data = movie_features)
ggqqplot(residuals(nv))

# Varianzhomogenität
leveneTest(data = movie_features, frame_position ~ hierarchy*duration)
plot(nv,1)


art_model <- art(frame_position ~ hierarchy * duration, data = movie_features)
m.art.anova = anova(art_model)

m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova

# mod_full <- lm(frame_position ~ hierarchy * duration, data = movie_features)
# mod_null <- lm(frame_position ~ 1, data = movie_features)  # intercept-only model
# anova(mod_null, mod_full)  # global F-test
