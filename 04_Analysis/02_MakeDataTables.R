# MAKE DATA TABLE FOR AFNI
# PREP ####
setwd('.../NIFTI')
data_dir = '.../NIFTI/'

subjects = c('BHM_437')

subject_nums = c(1)

subjects_sort = mark::sort_by(subjects, subject_nums)

Hierarchy = c('Shot', 'Scene')
num_hier = length(Hierarchy)
Duration = c('4s','12s','36s')
num_dur = length(Duration)

number_subejcts = length(subjects)
perCond = (number_subejcts * (number_subejcts -1))/2 # (N * (N-1))/2
pairs <- t(combn(number_subejcts, 2))  # 1128 pairs

# Generate all combinations of conditions
condition_grid <- expand.grid(Hierarchy = Hierarchy,
                              Duration = Duration,
                              KEEP.OUT.ATTRS = FALSE,
                              stringsAsFactors = FALSE)

# Repeat the subject pairs across conditions
subject_pairs <- data.frame(pairs[rep(1:nrow(pairs), each = num_hier * num_dur), ])
names(subject_pairs) = c("SubNum1", "SubNum2")

subject_pairs$Subj1 = subjects_sort[subject_pairs$SubNum1]
subject_pairs$Subj2 = subjects_sort[subject_pairs$SubNum2]

# Repeat the condition grid to match the number of pairs
conditions <- condition_grid[rep(1:(num_hier*num_dur), times = nrow(pairs)), ]

# Build the empty data frame
isc_table <- data.frame(
  Subj1 = subject_pairs$Subj1,
  Subj2 = subject_pairs$Subj2,
  Hierarchy = conditions$Hierarchy,
  Duration = conditions$Duration
)

curDate = format(Sys.Date(), "_%Y_%m_%d")

# Check the first few rows
head(isc_table)
isc_table$InputFile = NA

# create file names
for (i in 1:length(isc_table$Subj1)) {
  path = paste0(data_dir, isc_table$Subj1[i], '/ISC_Pair_maps_M13/ISC_Pair_FishZ_',
                isc_table$Hierarchy[i],'_',isc_table$Duration[i],'_',
                isc_table$Subj1[i],'_',isc_table$Subj2[i],'.nii \\')
  isc_table$InputFile[i] = path
}

# STEP 1 Condition-wise files ####
for (i in c('Shot', 'Scene')) {
  for (j in c('4s','12s','36s')){
    isc_out_tmp = isc_table[isc_table$Hierarchy == i & isc_table$Duration == j, c(1,2,3,4,5)]
    write.table(isc_out_tmp, file = paste0(data_dir,"isc_con_table_",i ,"_",j,curDate,".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
}

# STEP 2 Contrast maps ####
filesep = .Platform$file.sep
afniFolder = '_AFNI_Analysis'

# Contrasts
# 1.Shot vs Scene
# 2. 04s vs 12s
# 3. 04s vs 36s
# 4. 12s vs 36s

# all other possible maps, partially necessary for interactions
# 5. Shot04s vs Scene04s
# 6. Shot12s vs Scene12s
# 7. Shot36s vs Scene36s
# 8. Shot04s vs Shot12s
# 9. Shot04s vs Shot36s
# 10. Shot12s vs Shot36s
# 11. Scene04s vs Scene12s
# 12. Scene04s vs Scene36s
# 13. Scene12s vs Scene36s

contrasts = c( "01_Shot_Scene",
               "02_04s_12s",
               "03_04s_36s",
               "04_12s_36s",
               "05_Shot04_Scene04",
               "06_Shot12s_Scene12s",
               "07_Shot36s_Scene36s",
               "08_Shot04s_Shot12s",
               "09_Shot04_Shot36s",
               "10_Shot12s_Shot36s",
               "11_Scene04s_Scene12s",
               "12_Scene04s_Scene36s",
               "13_Scene12s_Scene36s")

outRoot = file.path(data_dir, afniFolder)

for (con in contrasts) {
  
  outDir = file.path(outRoot,con)
  if (!dir.exists(outDir)) {dir.create(outDir)}
  
  # Build the empty data frame
  isc_table <- data.frame(
    Subj1 = subject_pairs$Subj1,
    Subj2 = subject_pairs$Subj2,
    Contrast = con)
  
  isc_table$InputFile = NA
  
  for (i in 1:length(isc_table$Subj1)) {
    path = paste0(data_dir, filesep, isc_table$Subj1[i],filesep,con,filesep,
                  isc_table$Subj1[i],'_',isc_table$Subj2[i],'+orig \\')
    isc_table$InputFile[i] = path
  }
  
  write.table(isc_table, file = paste0(outDir,filesep,"isc_datatable_",con,".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# STEP 3 INTERACTION Contrast maps ####
# Different effect of Scene vs. Shot in different levels of Duration
# 1.(Shot04s vs Scene04s) vs (Shot12s vs Scene12s)
# 2.(Shot04s vs Scene04s) vs (Shot36s vs Scene36s)
# 3.(Shot12s vs Scene12s) vs (Shot36s vs Scene36s)


contrasts = c("INT01_12sSceneShot_04sSceneShot",
              "INT02_36sSceneShot_04sSceneShot",
              "INT03_36sSceneShot_12sSceneShot")

outRoot = file.path(data_dir, afniFolder)

for (con in contrasts) {
  
  outDir = file.path(outRoot,con)
  if (!dir.exists(outDir)) {dir.create(outDir)}
  
  # Build the empty data frame
  isc_table <- data.frame(
    Subj1 = subject_pairs$Subj1,
    Subj2 = subject_pairs$Subj2,
    Contrast = con)
  
  isc_table$InputFile = NA
  
  for (i in 1:length(isc_table$Subj1)) {
    path = paste0(data_dir, filesep, isc_table$Subj1[i],filesep,con,filesep,
                  isc_table$Subj1[i],'_',isc_table$Subj2[i],'+orig \\')
    isc_table$InputFile[i] = path
  }
  
  write.table(isc_table, file = paste0(outDir,filesep,"isc_datatable_",con,".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
