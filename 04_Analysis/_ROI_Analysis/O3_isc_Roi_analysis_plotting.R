# ISC ROI analysis ####
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(ggplot2)


study_dir = '.../NIFTI'
setwd(study_dir)

subjects = c("SHM_134","BST_790","KLN_831","SNE_929","UBN_644","MWL_955",
             "UHE_774","APN_733","SRE_111","CDT_864","COK_281","MDD_513",
             "MGH_898","SHG_243","CEN_612","DKL_934","KTF_128","SUA_130",
             "SHN_281","APN_759","MBZ_294","CHR_114","SCE_360","SRN_493",
             "NPN_211","NTA_213","DHR_253","GIN_112","JTR_424","MKN_944",
             "KWL_123","LMS_823","RRE_130","SBD_234","SIN_255","SLG_581",
             "RGH_912","BHG_483","AWR_117","KMR_411","BGU_132","CHM_519",
             "HMN_829","ECD_152","KBD_286","IWL_121","HDN_119","GCD_923")

subject_nums = c( 1, 2, 3, 6, 8, 9,
                 10,20,25,14,42,13,
                 15,11,35,17,45,34,
                  5,16,18,19,33,44,
                 38,37,21,22,23,24,
                 41,40,28,12,31,43,
                 30,27, 4,48,36,46,
                  7,32,39,47,26,29)

subjects_sort = mark::sort_by(subjects, subject_nums)

number_subjects = length(subjects_sort)

hier = c("Scene", "Shot")
dur = c("4s", "12s", "36s")

all_isc_list <- list()

for (i in 1:number_subjects) {
  print(subjects_sort[i])
  subj_dir = file.path(study_dir, subjects_sort[i])
  
    for (h in hier) {
      print(h)
      for (d in dur) {
        print(d)
          iscfile = paste0(subjects_sort[i],"_ISCpair_ROIs_", h, "_", d, ".csv")
          iscpath = file.path(subj_dir, iscfile)
          
          isc_tmp = read.csv(iscpath)
          
          # Pivot to long format
          isc_long = isc_tmp %>%
            pivot_longer(
              cols = names(isc_tmp)[5:52],            
              names_to = "Subj2",
              values_to = "ISC_pair_FishZ"
            ) %>%
            filter(!is.na(ISC_pair_FishZ))   
          
          isc_long$Hemisphere = ifelse(startsWith(isc_long$ROI, "Left"), "L", "R")
          isc_long$ROI_b = stringr::str_split_fixed(isc_long$ROI, "_", 2)[,2]
          
          all_isc_list[[length(all_isc_list) + 1]] = isc_long
      }
    }
}

iscData = bind_rows(all_isc_list)
summary(iscData)
iscData$Subj1 = factor(iscData$Subj1, levels = subjects_sort, labels = subjects_sort)
iscData$Subj2 = factor(iscData$Subj2, levels = subjects_sort, labels = subjects_sort)
iscData$Hierarchy = factor(iscData$Hierarchy, levels = c("Shot", "Scene"), labels = c("Shot", "Scene"))
iscData$Duration = factor(iscData$Duration, levels = dur, labels = dur)

iscData = read.csv(".../Shared_data/ROI_data.csv")

round(iscData$ISC_pair_FishZ[iscData$Subj1 == subjects_sort[1] & iscData$Subj2 == subjects_sort[2]],2) == 
round(iscData$ISC_pair_FishZ[iscData$Subj2 == subjects_sort[1] & iscData$Subj1 == subjects_sort[2]],2)

# Create a look up table for subject rank
subject_rank <- setNames(seq_along(subjects_sort), subjects_sort)

# Add rank columns to ISC data
iscData$rank1 <- subject_rank[iscData$Subj1]
iscData$rank2 <- subject_rank[iscData$Subj2]

# Keep only rows where Subj1 has a lower rank than Subj2
iscData_unique_brms <- iscData %>%
  filter(rank1 < rank2) %>%
  select(-rank1, -rank2)  # Clean up helper columns

iscData$rank1 = iscData$rank2 = NULL

save(iscData_unique_brms, file = "iscData_unique_brms.RData")

##LME4 models

# I need to use the fully crossed data, to accurately model the covariance structrue of the data
# To circumvent the impact of partially crossed random effect, here
# we loosen the index constraint of ij> in (2) to ij≠ , thereby
# incorporating both the upper and lower triangular portions of the Z
# ISC data (with redundancy) as input (Cheng et al., 2017)

####  Separate models per ROI 
all_rois = unique(iscData$ROI)

results_list = list()
for (roi in all_rois) {
  data_roi_sub = subset(iscData, iscData$ROI == roi)
  
  model = lmer(ISC_pair_FishZ ~ Hierarchy*Duration + (1|Subj1) + (1|Subj2),
               data = data_roi_sub)
  
  results_list[[roi]] = model
  
}


outDir = ".../ROI_analysis/_ANALYSIS/model_check_performance"

for (roi in all_rois) {
  png(file = file.path(outDir, paste0("check_model_",roi,".png")), width = 1200, height = 800)
  performPlot = performance::check_model(results_list[[roi]])
  print(performPlot)
  dev.off()
}


#####'STEP 2: Define a bootstrap function ######

#This will be passed to bootMer() and will return estimated marginal means for your fixed effects
#(e.g., for combinations of Hierarchy and Duration).
get_emm_estimates <- function(model) {
  emm = emmeans(model, ~ Hierarchy * Duration)
  output1 = as.numeric(summary(emm)$emmean)  # Only return the EMMs
  cont = contrast(emm, method = "pairwise")
  output2 = as.numeric(summary(cont)$estimate)  # returns vector of contrast estimate
  output = c(output1, output2)
  return(output)
}

emm1 =  emmeans(results_list$Left_Angular_Gyrus_Anatomy_PG, ~ Hierarchy * Duration)
emm1_sum = summary(emm1)
contrast1 = summary(contrast(emm1, method = "pairwise"))
contrasts_con1 = contrast1$contrast
# Replace " - " with "_"
con_names <- gsub(" - ", "_", contrasts_con1)
# Replace spaces " " with "_"
con_names <- gsub(" ", "_", con_names)

names_emm = paste0(emm1_sum$Hierarchy,'_', emm1_sum$Duration) 

all_names_boot = c(names_emm, con_names)

get_emm_estimates(results_list$Left_Angular_Gyrus_Anatomy_PG)

library(pbapply)
set.seed(123)  # for reproducibility

boot_results = list()

boot_reps = 500  # You can increase this if you have time

boot_results = pblapply(names(results_list), function(roi) {
  
  model = results_list[[roi]]
  
  # bootMer with emmeans-based function
  boot = bootMer(
    model,
    FUN = get_emm_estimates,
    nsim = boot_reps,
    use.u = TRUE,
    type = "parametric",
    parallel = "multicore",  # use "snow" on Windows
    ncpus = 4
  )
  
  list(
    ROI = roi,
    model = model,
    boot = boot
  )
})

save.image(".../isc_Roi_analysis_2.RData")

# Function to get summary stats per ROI
summarize_boot_emm <- function(boot, roi) {
  boot_matrix = boot$t[,1:6]
  colnames(boot_matrix) = names_emm
  
  # Get bootstrap means and 95% CIs
  apply(boot_matrix, 2, function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      lower = quantile(x, 0.025, na.rm = TRUE),
      upper = quantile(x, 0.975, na.rm = TRUE)
    )
  }) |> t() |> as.data.frame()
}

# Function to get Contrasts  stats per ROI
summarize_boot_con <- function(boot, roi) {
  boot_matrix = boot$t[,7:21]
  colnames(boot_matrix) = con_names
  
  # Get bootstrap means and 95% CIs
  output = apply(boot_matrix, 2, function(x) {
    c(mean = mean(x, na.rm = TRUE),
      lower = quantile(x, 0.025, na.rm = TRUE),
      upper = quantile(x, 0.975, na.rm = TRUE)
    )
  }) |> t() |> as.data.frame()
  
  # calcualte two-sided p-values
  output$p_value <- apply(
    boot_matrix,  # only contrast columns
    2,
    function(x) {
      2 * min(mean(x <= 0), mean(x >= 0))
    }
  )
  return(output)
}


# Combine all ROI summaries
boot_summary_emm = do.call(rbind, lapply(seq_along(boot_results), function(i) {
  roi = boot_results[[i]]$ROI
  summ = summarize_boot_emm(boot_results[[i]]$boot, roi)
  summ$ROI = roi
  summ$Condition = names_emm
  summ
}))

# Combine all contrast summaries
boot_summary_con = do.call(rbind, lapply(seq_along(boot_results), function(i) {
  roi = boot_results[[i]]$ROI
  summ = summarize_boot_con(boot_results[[i]]$boot, roi)
  summ$ROI = roi
  summ$Condtion = con_names
  summ
}))

save.image(".../isc_Roi_analysis_3.RData")


# Correct the p-Values
boot_summary_con$p_value_FDR = p.adjust(boot_summary_con$p_value, method = "BH")
boot_summary_con$p_value_FDRy = p.adjust(boot_summary_con$p_value, method = "BY")
alpha = 0.05
boot_summary_con$sig_FDR = ifelse(boot_summary_con$p_value_FDR < alpha, 1, 0)
boot_summary_con$sig_FDRy = ifelse(boot_summary_con$p_value_FDRy < alpha, 1, 0)

outDir = ".../Movie-HINTS_Data/ROI_analysis"

writexl::write_xlsx(boot_summary_con, path = file.path(outDir, "Boot_summary_all_Contrasts.xlsx"))

# Plot the comparisions with bootstrap results

rois_clean = c("Left Angular Gyrus",
               "Left Heschl's Gyrus (A1)",
               "Left Calcarine Sulcus (V1)",
               "Left PCC/Precuneus",
               "Left hV4",
               "Right Angular Gyrus",
               "Right Heschl's Gyrus (A1)",
               "Right Calcarine Sulcus (V1)",
               "Right PCC/Precuneus",
               "Right hV4",
               
               "Left Parahippocampal Gyrus",
               "Left TPJ",
               "Left anterior Hippocampus",
               "Left mPFC",
               "Left posterior Hippocampus",
               "Right Parahippocampal Gyrus",
               "Right TPJ",
               "Right anterior Hippocampus",
               "Right mPFC",
               "Right posterior Hippocampus")

rois_clean_abb = c("left AG",
                   "left A1",
                   "left V1",
                   "left pmDMN",
                   "left hV4",
                   "righ AG",
                   "right A1",
                   "right V1",
                   "right pmDMN",
                   "right hV4",
                   
                   "left PHG",
                   "left TPJ",
                   "left aHC",
                   "left mPFC",
                   "left pHC",
                   "right PHG",
                   "right TPJ",
                   "right aHC",
                   "right mPFC",
                   "right pHC")


colors = list(c("#978802","#fde935"),
              c("#2d8630","#66cc69"),
              c("#135350","#26a6a0"),
              c("#1f2a47","#3d548f"),
              c("#280132","#500264"))

colorPick = c(4,1,2,5,3,4,1,2,5,3,
              1,4,3,5,2,1,4,3,5,2)

names(boot_summary_emm)[2:3] = c("lower_0025", "upper_0975")

conditions = str_split_fixed(boot_summary_emm$Condition,'_', 2)
boot_summary_emm$Hierarchy = conditions[,1]
boot_summary_emm$Duration = conditions[,2]

boot_summary_emm$Hierarchy = factor(boot_summary_emm$Hierarchy, levels = c("Shot", "Scene"), labels = c("Shot", "Scene"))
boot_summary_emm$Duration = factor(boot_summary_emm$Duration, levels = c("4s", "12s", "36s"), labels = c("4s", "12s", "36s"))

outDir = ".../ROI_analysis/_ANALYSIS/boot_plots/_color_picked_no_legend"
text_size= 8
for (roi in 1:length(all_rois)) {
  boot_df = data.frame(subset(boot_summary_emm, boot_summary_emm$ROI == all_rois[roi]))
  
  gg = ggplot(boot_df, aes(x = Duration, y = mean, color = Hierarchy, group = Hierarchy)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(aes(ymin = lower_0025, ymax = upper_0975), width = 0.08, linewidth = 0.8) +
    theme_classic() + ggtitle(rois_clean_abb[roi]) + ylab("Average ISC") + scale_color_manual(values = colors[[colorPick[roi]]])+
    labs(color = "Content Level") + theme(legend.position="none") +
    theme(text = element_text(size = text_size, colour ="black"), axis.text.x=element_text(size = text_size, colour="black"), 
          axis.text.y=element_text(size = text_size, colour="black"))
  
  ggsave(filename = paste0("boot_", all_rois[roi], ".png"),
         plot = gg,  path = outDir, width = 800, height = 800, units = "px", dpi = 300)
  
}
