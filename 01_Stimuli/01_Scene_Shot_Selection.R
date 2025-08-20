# DEFINTIONS & LIBRARIES
setwd(".../studyforrest/gump_emotions")

library(dplyr)
library(stringr)
library(ggplot2)
library(writexl)

scenes = read.csv("movie_scenes.csv", header = F)
shots = read.csv("movie_shots.csv", header = F)
names(scenes) = c("seconds_raw", "place", "time", "location")
names(shots) = c("seconds_raw")

scene_markers = read.csv(".../AdditionalFiles/scene_timestamps_20241031.csv")
shot_markers = read.csv(".../AdditionalFiles/shot_timestamps_20241112.csv")

emotions = read.delim(file = "./segmentation/emotions_av_shots_thr50.tsv", header = FALSE)
names(emotions) = c("startTime", "endTime", "desc")

# REMOVE DUPLICATES
filter4duplicates = matrix(TRUE, length(shots$seconds_raw),1)
for (i in 2:length(shots$seconds_raw)) {
  if (shots$seconds_raw[i-1] == shots$seconds_raw[i]) {
    filter4duplicates[i] = FALSE
  }
}
shots = data.frame(shots$seconds_raw[filter4duplicates])
names(shots) = c("seconds_raw")


# ADD START AND END TIME IN CORRECT FORMAT
scenes$hours = floor(scenes$seconds_raw/3600)
scenes$minutes = floor(scenes$seconds_raw %% 3600 /60)
scenes$seconds = scenes$seconds_raw %% 60

shots$hours = floor(shots$seconds_raw/3600)
shots$minutes = floor(shots$seconds_raw %% 3600 /60)
shots$seconds = shots$seconds_raw %% 60

library(stringr)
ending = 7082.20
floor(ending/3600)
floor(ending %% 3600 /60)
ending %% 60

scenes$start_time = paste0(str_pad(scenes$hours, 2, pad = "0"), ':', str_pad(scenes$minutes, 2, pad = "0"), ':', sprintf(scenes$seconds, fmt = '%05.2f'))

shots$start_time = paste0(str_pad(shots$hours, 2, pad = "0"), ':', str_pad(shots$minutes, 2, pad = "0"), ':', sprintf(shots$seconds, fmt = '%05.2f'))

scenes$end_time = c(scenes$start_time[2:length(scenes$start_time)],"01:58:02.20")
shots$end_time = c(shots$start_time[2:length(shots$start_time)],"01:58:02.20")


# SCENE DURATION AND SHOT NUMBER
shots$duration = c(shots$seconds_raw[2:length(shots$seconds_raw)], ending) - shots$seconds_raw
scenes$duration = c(scenes$seconds_raw[2:length(scenes$seconds_raw)],ending)  - scenes$seconds_raw

scenes$is_shot_start_raw = scenes$seconds_raw %in% shots$seconds_raw

scenes$seconds_Round = round(scenes$seconds_raw,0)
shots$seconds_Round = round(shots$seconds_raw,0)

scenes$is_shot_start_Round = scenes$seconds_Round %in% shots$seconds_Round

scenes$howManyShots = NA
shots$is_seconds_Round = shots$seconds_Round %in% scenes$seconds_Round
scenes$seconds_Round_corrected = scenes$seconds_Round

for (i in 1:(length(scenes$seconds_raw)-1)) {
  if (scenes$is_shot_start_Round[i] == TRUE){
    scenes$howManyShots[i] = max(which(shots$seconds_Round <= scenes$seconds_Round[i+1])) -  which(shots$seconds_Round == scenes$seconds_Round[i])
  }
}

## manually correct fro WW2 etc.
scenes$seconds_Round_corrected[scenes$is_shot_start_Round == FALSE | scenes$howManyShots == 0 | is.na(scenes$howManyShots)] =
  c(0, 1136, 1145, 1856, 2312, 2313, 5474, 5482, 6945)
# old:   17 1136 1144 1857 2312 2313 5474 5481 6945
# WorldWar2 2313, ist als shot vergessen
# letzte Szene hat 19 shots
shots$is_seconds_Round_corrected = shots$seconds_Round %in% scenes$seconds_Round_corrected
scenes$is_shot_start_Round_corrected = scenes$seconds_Round_corrected %in% shots$seconds_Round
scenes$howManyShots_corrected = NA


for (i in 1:(length(scenes$seconds_raw)-1)) {
  if (scenes$is_shot_start_Round_corrected[i] == TRUE){
    scenes$howManyShots_corrected[i] = max(which(shots$seconds_Round <= scenes$seconds_Round_corrected[i+1])) -  which(shots$seconds_Round == scenes$seconds_Round_corrected[i])
  }
}

# last manual touchup: correcting WW1 and WW2 and last scene
scenes$howManyShots_corrected[scenes$howManyShots_corrected == 0 | is.na(scenes$howManyShots_corrected)] = c(1,1,19)

# SCENE SELECTION ####
scenes4select = data.frame(start_time = scene_markers$start_time,
                           place = scenes$place,
                           time = scenes$time,
                           location = scenes$location, 
                           duration = scenes$duration, 
                           howManyShots = scenes$howManyShots_corrected, 
                           stringsAsFactors = F)

# CONNECT SHOTS TO SCENE
shots$whichScene = 0
m = 0
for (i in 1:length(shots$whichScene)) {
  if (shots$is_seconds_Round_corrected[i]) {
    m = m+1
  }
  shots$whichScene[i] = m
}
# correct for WW2 scene
shots$whichScene[293:length(shots$whichScene)] = shots$whichScene[293:length(shots$whichScene)] +1
shots$seconds_round1 = round(shots$seconds_raw, 1)


# ADD EMOTIONAL AROUSAL
emotions$startTime %in% shots$seconds_round1
emotions$endTime %in% emotions$startTime

emotions <- emotions %>%
  mutate(
    person = str_extract(emotions$desc, "(?<=char=)[^\\s]+"),
    arousal = as.numeric(str_extract(emotions$desc, "(?<=arousal=)-?\\d+\\.\\d+")),
    val_pos = as.numeric(str_extract(emotions$desc, "(?<=val_pos=)-?\\d+\\.\\d+")),
    val_neg = as.numeric(str_extract(emotions$desc, "(?<=val_neg=)-?\\d+\\.\\d+"))
  )

# FORRESTVO ist not a person
emotions$person[emotions$person == "FORRESTVO"] = NA

emotionsAgg = aggregate(cbind(arousal, val_pos, val_neg) ~ startTime, data = emotions, mean)

emotionsAgg$peoplePerShot = 0
people= list()

for (i in unique(emotions$startTime)) {
  
  peopleInShot = list((emotions$person[emotions$startTime == i]))
  people = append(people, peopleInShot)
  
  emotionsAgg$peoplePerShot[emotionsAgg$startTime == i] = length(peopleInShot[[1]])
}

emotionsAgg$people = people

shots$arousal = NA
shots$val_pos = NA
shots$val_neg = NA
shots$peoplePerShot = NA
shots$people = NA
for (i in 1:length(emotionsAgg$startTime)) {
  shots$arousal[shots$seconds_round1 == emotionsAgg$startTime[i]] = emotionsAgg$arousal[i]
  shots$val_pos[shots$seconds_round1 == emotionsAgg$startTime[i]] = emotionsAgg$val_pos[i]
  shots$val_neg[shots$seconds_round1 == emotionsAgg$startTime[i]] = emotionsAgg$val_neg[i]
  shots$peoplePerShot[shots$seconds_round1 == emotionsAgg$startTime[i]] = emotionsAgg$peoplePerShot[i]
  shots$people[shots$seconds_round1 == emotionsAgg$startTime[i]] = emotionsAgg$people[i]
}

sum(shots$seconds_round1 %in% emotionsAgg$startTime)
sort(table(shots$seconds_round1[shots$seconds_round1 %in% emotionsAgg$startTime]))


emoScene = aggregate(cbind(arousal,val_pos,val_neg) ~ whichScene, data = shots, mean)
emoScene = cbind(emoScene, aggregate(cbind(arousal,val_pos,val_neg) ~ whichScene, data = shots, sd))
emoScene[,5] = NULL
names(emoScene)[2:7] = c("arousal_mean", "val_pos_mean", "val_neg_mean", "arousal_sd", "val_pos_sd", "val_neg_sd")

# ADD HOW MANY PEOPLE ARE IN A SCENE
emoScene$peoplePerScene = 0
peopleScene= list()

for (i in unique(emoScene$whichScene)) {
  
  peopleInScene = shots$people[shots$whichScene == i]
  peopleInScene = unlist(peopleInScene)
  peopleInScene = unique(peopleInScene[!is.na(peopleInScene)])
  peopleScene = append(peopleScene, list(peopleInScene))
  
  emoScene$peoplePerScene[emoScene$whichScene == i] = length(peopleInScene)
}

emoScene$people = peopleScene

scenes4select$peoplePerScene = NA
scenes4select$arousal_mean = NA
scenes4select$val_pos_mean = NA
scenes4select$val_neg_mean = NA
scenes4select$arousal_sd = NA
scenes4select$val_pos_sd = NA
scenes4select$val_neg_sd = NA

scenes4select$peoplePerScene[emoScene$whichScene] = emoScene$peoplePerScene
scenes4select$peoplePerScene[scenes4select$peoplePerScene == 0] = NA
scenes4select$arousal_mean[emoScene$whichScene] = emoScene$arousal_mean
scenes4select$val_pos_mean[emoScene$whichScene] = emoScene$val_pos_mean
scenes4select$val_neg_mean[emoScene$whichScene] = emoScene$val_neg_mean
scenes4select$arousal_sd[emoScene$whichScene] = emoScene$arousal_sd
scenes4select$val_pos_sd[emoScene$whichScene] = emoScene$val_pos_sd
scenes4select$val_neg_sd[emoScene$whichScene] = emoScene$val_neg_sd



## SHOT FOCUS #####

# 1. Shots with People in them 
shots2 = shots
shots2$description = shot_markers$Description 

split_text = str_split(shots2$description, " ")

# Extract the last three words and concatenate the rest
result <- lapply(split_text, function(words) {
  n <- length(words)
  if (n >= 3) {
    last_three <- words[(n - 2):n]  # Get last three words
    rest <- paste(words[1:(n - 3)], collapse = " ")  # Concatenate the rest
  } else {
    # If fewer than 3 words
    last_three <- c(rep(NA, 3 - n), words)  # Pad with NAs to make 3
    rest <- ""  # No "rest"
  }
  list(last_three = last_three, rest = rest)
})

# Convert the result into columns
shots2$place <- sapply(result, function(x) x$rest)
shots2$time <- sapply(result, function(x) x$last_three[1])
shots2$location <- sapply(result, function(x) x$last_three[2])
shots2$repetition_location <- sapply(result, function(x) x$last_three[3])

shots4select = data.frame(start_time = shot_markers$start_time,
                          place = shots2$place,
                          time = shots2$time,
                          location = shots2$location,
                          duration = shots2$duration,
                          howManyShots = 1,
                          arousal_mean = shots2$arousal,
                          val_pos_mean = shots2$val_pos,
                          val_neg_mean = shots2$val_neg, 
                          arousal_sd = NA,
                          val_pos_sd = NA,
                          val_neg_sd = NA,
                          peoplePerScene = shots2$peoplePerShot,
                          people = shots2$repetition_location)

shots4select$time[shots4select$time == "DAWN"] = "DAY"
shots4select$location[shots4select$location == "DAY"] = "EXT"

shotSelect1 = shots4select[!(is.na(shots4select$peoplePerScene)),]


# more manual adjustment von shots und scenes
writexl::write_xlsx(shots4select, "shots4select.xlsx")
writexl::write_xlsx(scenes4select, "scenes4select.xlsx")
cumSumshots = data.frame(c(1, cumsum(scenes4select$howManyShots[1:204])+1))
writexl::write_xlsx(cumSumshots, "cumSumShots.xlsx")
cumSumshots[68:205,1] = cumSumshots[68:205,1] -1


shots4select2 = readxl::read_xlsx(".../AdditionalFiles/shots4select.xlsx")
shots4select2$`0` = NULL
shots4select2$whichScene = shots$whichScene
### die label waren nicht richtig zugeordnet bei den Shots


# Function to convert HH:MM:SS:FF to total frames
time_to_frames <- function(time_str, fps = 25) {
  # Split the time string into HH, MM, SS, FF
  parts <- strsplit(time_str, ":")[[1]]
  hh <- as.numeric(parts[1])
  mm <- as.numeric(parts[2])
  ss <- as.numeric(parts[3])
  ff <- as.numeric(parts[4])
  
  # Convert to total frames
  total_frames <- (hh * 3600 + mm * 60 + ss) * fps + ff
  return(total_frames)
}

# Calculate frame differences
frame_differences <- mapply(function(t1, t2) {
  abs(time_to_frames(t1) - time_to_frames(t2))
}, scenes4select$start_time, shots4select2$start_time[cumSumshots[,1]])

# Print the differences
frame_differences # maximally 4 frame differnces --> shots in cut version so we use thems


scenes4select$start_time_shots = shots4select2$start_time[cumSumshots[,1]]
scenes4select$start_time_shots[68] = scenes4select$start_time[68] # correct WW2
scenes4select$end_time_shots = c(scenes4select$start_time_shots[2:length(scenes4select$start_time)], scenes$end_time[length(scenes4select$start_time)])

# get the first or last shot of a scene
start_last_shot_of_scene = shots4select2$start_time[which(shots4select2$start_time %in% scenes4select$end_time_shots) -1]
start_last_shot_of_scene = c(start_last_shot_of_scene[1:67],scenes4select$start_time_shots[68] ,start_last_shot_of_scene[68:length(start_last_shot_of_scene)], shots4select2$start_time[length(shots4select2$start_time)])
scenes4select$start_last_shot_of_scene =start_last_shot_of_scene

frame_differences <- mapply(function(t1, t2) {
  abs(time_to_frames(t1) - time_to_frames(t2))
}, scenes4select$start_time, scenes4select$start_time_shots)
frame_differences


# get all shots2 that are
## 1. showing a person
## 2. a single shot scene
## 3. the first or last shot of a scene with at least 4 shots

# have people
shotSelect1_1 = shots4select2[!(is.na(shots4select2$peoplePerScene)),]

# are part of one shot scenes
shotSelect1_2 = shotSelect1_1[shotSelect1_1$whichScene %in% as.numeric(rownames(scenes4select[(!is.na(scenes4select$peoplePerScene)) & 
                                                                                                scenes4select$howManyShots == 1,])),]

# SCENE FOCUS
# 1. SELECTION SCENES WITH PEOPLE IN THEM that have at least 2 shots
scenesSelect1 = scenes4select[(!is.na(scenes4select$peoplePerScene)) & scenes4select$howManyShots > 1,]

summary(shotSelect1_2[c(3:7,13)])
summary(scenesSelect1[c(3:6,11,10)])

# Scenes have longer Duration, higehr arousal and more people
# remove extreme values from scenes
boxScene = boxplot(scenesSelect1$duration)
SceneDurationThresh = min(boxScene$out)
ShotDurationThresh = max(shotSelect1_2$duration)

boxScene = boxplot(scenesSelect1$arousal_mean)
# SceneArousalThresh = min(boxScene$out) # no outliers
ShotArousalThresh = max(shotSelect1_2$arousal_mean)

boxScene = boxplot(scenesSelect1$peoplePerScene)
# ScenePeopleThresh = min(boxScene$out) # no outliers
ShotPeopleThresh = max(shotSelect1_2$peoplePerScene)

scenesSelect2 = scenesSelect1[scenesSelect1$duration <= ShotDurationThresh & 
                                scenesSelect1$arousal_mean <=  ShotArousalThresh &
                                scenesSelect1$peoplePerScene <=  ShotPeopleThresh,]

whichScenesSelect2 = as.numeric(rownames(scenesSelect2))
shotSelect2 = shotSelect1_1[!(shotSelect1_1$whichScene %in% whichScenesSelect2),]

summary(shotSelect2[c(3:8,14)])
summary(scenesSelect2[c(2:7,13)])

SceneDurationThresh = 4.00
SceneArousalThresh = min(scenesSelect2$arousal_mean)
ScenePeopleThresh = min(scenesSelect2$peoplePerScene)

# only shots that have a similar min and max to the scene selection
shotSelect3 = shotSelect2[shotSelect2$duration >= SceneDurationThresh & 
                            shotSelect2$arousal_mean >=  SceneArousalThresh &
                            shotSelect2$duration <= ShotDurationThresh & 
                            shotSelect2$arousal_mean <=  ShotArousalThresh &
                            shotSelect2$peoplePerScene >=  ScenePeopleThresh,]

summary(shotSelect3[c(3:8,14)])
summary(scenesSelect2[c(2:7,13)])

# only shots that are part of 1 shot scenes
shotSelect4_1 = shotSelect3[shotSelect3$whichScene %in% as.numeric(rownames(scenes4select[(!is.na(scenes4select$peoplePerScene)) & 
                                                                                            scenes4select$howManyShots == 1,])),]


save.image(".../studyforrest/01_Scene_Shot_Selection_01.RData")

# Subsampling
scenesSelectSet = scenesSelect2[,c(2,3,4,5,11,10,6,7,8)]
shotsSelectSet = shotSelect4_1[,c(2,3,4,5,7,13,6,1,15)]

scenesSelectSet$Ratio4 = scenesSelectSet$duration/4
scenesSelectSet$Ratio12 = scenesSelectSet$duration/12
scenesSelectSet$Ratio36 = scenesSelectSet$duration/36

shotsSelectSet$Ratio4 = shotsSelectSet$duration/4
shotsSelectSet$Ratio12 = shotsSelectSet$duration/12
shotsSelectSet$Ratio36 = shotsSelectSet$duration/36

scenesSelectSet$whichScene = rownames(scenesSelectSet)

writexl::write_xlsx(scenesSelectSet, "scenesSelectSet.xlsx")
writexl::write_xlsx(shotsSelectSet, "shotsSelectSet.xlsx")


medAll = median(c(scenesSelectSet[,10], scenesSelectSet[,11], scenesSelectSet[,12], 
       shotsSelectSet$Ratio4, shotsSelectSet$Ratio12,shotsSelectSet$Ratio36))
box = boxplot(c(scenesSelectSet[,10], scenesSelectSet[,11], scenesSelectSet[,12], 
          shotsSelectSet$Ratio4, shotsSelectSet$Ratio12,shotsSelectSet$Ratio36))



library(clue)
## optimization of distance
# Match shots to scenes
match_scenes_to_shots <- function(shotsSelection, scenesSelection, vars, num_matches = nrow(shots), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  shots_1 = shotsSelection
  scenes_1 = scenesSelection
  
  # Compute pairwise distances between shots and scenes
  distance_matrix <- as.matrix(dist(rbind(shots_1[vars], scenes_1[vars])))
  shot_count <- nrow(shots_1)
  scene_count <- nrow(scenes_1)
  
  # Extract the distance matrix for shots vs. scenes
  distance_matrix <- distance_matrix[1:shot_count, (shot_count + 1):(shot_count + scene_count)]
  
  # Solve the assignment problem (minimizing total distance)
  match <- solve_LSAP(distance_matrix)
  
  # Create pairs
  pairs <- data.frame(
    shot_id = rownames(shotsSelection),
    scene_id = rownames(scenesSelection)[match],
    distance = distance_matrix[cbind(1:shot_count, match)]
  )
  
  # Add shot and scene details for matched pairs
  pairs <- cbind(
    pairs,
    shotsSelection[pairs$shot_id, vars, drop = FALSE],
    scenesSelection[pairs$scene_id, vars, drop = FALSE]
  )
  colnames(pairs) <- c("shot_id", "scene_id", "distance", 
                       paste0("shot_", vars), paste0("scene_", vars))
  
  # Identify remaining scenes
  remaining_scenes <- scenesSelection[-match, ]
  
  return(list(pairs = pairs, remaining_scenes = remaining_scenes, distance_matrix = distance_matrix))
}

# Example Usage
# Example datasets
shotsSelectSet_copy = shotsSelectSet
scenesSelectSet_copy = scenesSelectSet

# Match shots to scenes
vars <- c("duration", "arousal_mean", "peoplePerScene")
result <- match_scenes_to_shots(shotsSelectSet_copy, scenesSelectSet_copy, vars)

# Output
paired_scenes <- result$pairs
remaining_scenes <- result$remaining

# without normalization
result_nonNorm <- match_scenes_to_shots(shotsSelectSet_copy, scenesSelectSet_copy, vars)
paired_nonNorm_scenes <- result_nonNorm$pairs
remaining_nonNorm_scenes <- result_nonNorm$remaining
scenesCutinSelect = scenesSelect1[!(rownames(scenesSelect1) %in% rownames(scenesSelect2)),]

# all reamining possible shots from all scenes not matched yet
shotsSelect4_2 = shotSelect3[shotSelect3$whichScene %in% 
                             c(as.numeric(rownames(remaining_nonNorm_scenes)), 
                               as.numeric(rownames(scenesCutinSelect))),]


# find least likely to match scenes
least_likely_to_match <- function(remaining_scenes, shotsSelections, vars) {

  shotsSelections <- shotsSelections
  remaining_scenes <- remaining_scenes 
  
  # Compute pairwise distances between remaining scenes and shots
  distance_matrix <- outer(
    1:nrow(remaining_scenes), 
    1:nrow(shotsSelections), 
    Vectorize(function(i, j) {
      sum((remaining_scenes[i, vars] - shotsSelections[j, vars])^2)
    })
  )
  
  # Compute metrics for "hard-to-match" scenes
  avg_distances <- rowMeans(distance_matrix)
  min_distances <- apply(distance_matrix, 1, min)
  
  # Create a summary data frame
  results <- data.frame(
    scene_id = rownames(remaining_scenes),
    avg_distance = avg_distances,
    min_distance = min_distances
  )
  
  # Rank scenes by average distance (higher = harder to match)
  results <- results[order(-results$avg_distance), ]
  
  return(results)
}

remaining_results <- least_likely_to_match(remaining_nonNorm_scenes, shotsSelect4_2, vars)
print(head(remaining_results))

shotsSelect4_2_Set = shotSelect3[shotSelect3$whichScene %in% 
                               c(as.numeric(remaining_results$scene_id[1:2]), 
                                 as.numeric(rownames(scenesCutinSelect))),]

scenesSelectRemain = remaining_nonNorm_scenes[!(rownames(remaining_nonNorm_scenes) %in% remaining_results$scene_id[1:2]),]



remaining_scenes = scenesSelectRemain
secondary_shots = shotsSelect4_2_Set
which_scene_col = "whichScene"

#match_shots_to_remaining_scenes <- function(remaining_scenes, secondary_shots, vars, which_scene_col, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  remaining_scenes_i = remaining_scenes
  secondary_shots_i = secondary_shots
  scene_id_i = rownames(remaining_scenes_i)
  shot_id_i = rownames(secondary_shots_i)
  pairs_final = c()
  
  while(!is.null(remaining_scenes_i)) {
      

  # Compute the distance matrix
  distance_matrix <- as.matrix(dist(rbind(remaining_scenes_i[vars], secondary_shots_i[vars])))
  scenes_count <- nrow(remaining_scenes_i)
  shot_count <- nrow(secondary_shots_i)
  
  # Extract the distance matrix for shots vs. scenes
  distance_matrix <- distance_matrix[1:scenes_count, (scenes_count + 1):(scenes_count + shot_count)]
  
  
  # Solve the assignment problem (minimizing total distance)
  match <- solve_LSAP(distance_matrix)# Solve the assignment problem
  
  
  # Create pairs
  pairs <- data.frame(
    scene_id = scene_id_i,
    shot_id = shot_id_i[match],
    distance = distance_matrix[cbind(1:scenes_count, match)]
  )
  
  pairs$whichScene = unlist(secondary_shots[pairs$shot_id, which_scene_col])
  unselect = c()
  for (j in unique(pairs$whichScene)) {
    pairsub = pairs[pairs$whichScene == j,]
    if(dim(pairsub)[1] > 1){
      pairsub = pairsub[order(pairsub$distance),]
      unselect = c(unselect, rownames(pairsub)[2:(dim(pairsub)[1])])
    }
  }
  
  if (!is.null(unselect)) {
    pairsSelect = pairs[-as.numeric(unselect),]
  pairs_final = rbind(pairs_final, pairsSelect)
  remaining_scenes_i = remaining_scenes_i[as.numeric(unselect),]
  scene_id_i = scene_id_i[as.numeric(unselect)]
  shot_id_i = shot_id_i[!(secondary_shots_i$whichScene %in% pairsSelect$whichScene)]
  secondary_shots_i = secondary_shots_i[!(secondary_shots_i$whichScene %in% pairsSelect$whichScene), ]
  } else {
    pairs_final = rbind(pairs_final, pairs)
    remaining_scenes_i = NULL; scene_id_i = NULL
    shot_id_i = shot_id_i[!(secondary_shots_i$whichScene %in% pairs$whichScene)]
    secondary_shots_i = secondary_shots_i[!(secondary_shots_i$whichScene %in% pairs$whichScene), ]
  }
  

  } 
#   return(pairsfinal)
# }

# results_remainding_scenes = match_shots_to_remaining_scenes(scenesSelectRemain, shotsSelect4_2_Set, vars, "whichScene")
  
  pairsFull <- cbind(
    pairs_final,
    shotsSelect4_2_Set[pairs_final$shot_id, vars, drop = FALSE],
    scenesSelectRemain[pairs_final$scene_id, vars, drop = FALSE]
  )
  colnames(pairsFull) <- c("scene_id", "shot_id", "distance", "whichScene",
                       paste0("shot_", vars), paste0("scene_", vars))

length(unique(c(paired_nonNorm_scenes$scene_id, shotsSelectSet_copy$whichScene, pairsFull$scene_id, pairsFull$whichScene)))

scenes_paired_id = c(paired_nonNorm_scenes$scene_id, pairsFull$scene_id)
scenesSelectSet_paired = scenesSelectSet_copy[scenes_paired_id,]

shotsSelectSet_paried = rbind(shotsSelectSet_copy[paired_nonNorm_scenes$shot_id,1:9], 
                              shotsSelect4_2_Set[pairsFull$shot_id,c(2,3,4,5,7,13,6,1,15)])
save(scenesSelectSet_paired, shotsSelectSet_paried, file = "SelectSet_paired_20241130.RData" )  


### scene selection
# Define  function to split the dataset



split_dataset_balanced <- function(data, sizes, variables, global_means, global_sds, 
                                   factors, max_repeats = 1000, propThresh = 0.1) {
  subsets <- list()  # To store subsets
  
  # Compute global proportions of factor levels
  facotors_glob = global_props = list()
  for (f in seq_along(factors)) {
    facotors_glob[[factors[f]]] = factor(data[[factors[f]]], levels = unique(data[[factors[f]]]))
    global_props[[factors[f]]] = prop.table(table(facotors_glob[[factors[f]]]))
  }
  
  #print(global_props)
  
  for (size in sizes) {
    repeat {
      # Randomly sample rows
      sampled_indices <- sample(1:nrow(data), size)
      subset <- data[sampled_indices, ]
      
      # Check if the subset's means are within 1 SD of the global means
      subset_means <- colMeans(subset[, variables])
      mean_diff <- abs(subset_means - global_means)
      #print(mean_diff <= global_sds) # for finding issues
      within_one_sd <- all(mean_diff <= global_sds)
      
      
      # Check if the factor distributions are balanced
      factor_balanced <- TRUE
      for (i in seq_along(factors)) {
        factor <- factors[i]
        #print(factor)
        #print(factor(subset[[factor]], levels = levels(facotors_glob[[factor]])))
        subset_props <- prop.table(table(factor(subset[[factor]], levels = levels(facotors_glob[[factor]]))))
       #  print(subset_props)
        prop_diff <- abs(subset_props - global_props[[i]])

      #print(prop_diff) # for finding issues in convergence
        if (any(prop_diff > propThresh)) {  # Tolerance for imbalance
          factor_balanced <- FALSE
          break
        }
      }
      
      # If both conditions are met, accept the subset
      if (within_one_sd && factor_balanced) {
        subsets[[length(subsets) + 1]] <- subset
        data <- data[-sampled_indices, ]  # Remove selected rows from data
        break
      }
      
      # Decrement the repeat counter
      max_repeats <- max_repeats - 1
      if (max_repeats <= 0) {
        stop("Could not find a suitable subset within the maximum number of repeats.")
      }
    }
  }
  
  return(subsets)
}

# Step 2: Compute the global means and SDs
scene_global_means <- colMeans(scenesSelectSet_paired[, c("duration", "arousal_mean", "peoplePerScene", "howManyShots")])
scene_global_sds <- apply(scenesSelectSet_paired[, c("duration", "arousal_mean", "peoplePerScene", "howManyShots")], 2, sd)
scene_factors = c("time", "location")

# Step 3: Split the data into subsets
framesBlock = 648

sizes <- c(framesBlock/108, framesBlock/36, framesBlock/12)  # Desired sizes
vars <- c("duration", "arousal_mean", "peoplePerScene", "howManyShots")  # Variables to balance

# global_sd are too generous hard code preferred range
scene_global_range = scene_global_sds/c(4,8,8,8)

# select shots
subsetsScenes <- split_dataset_balanced(scenesSelectSet_paired, sizes, vars, scene_global_means, scene_global_sds, scene_factors, max_repeats = 10000, propThresh = 0.1)

data = scenesSelectSet_paired
global_means = scene_global_means
global_sds = scene_global_sds
factors = scene_factors
max_repeats = 100
propThresh = 0.1

optimize_balanced_split <- function(data, sizes, vars, global_means, global_sds, factors, max_repeats = 1000, step = 0.1) {
  # Initialize the ranges
  current_sds <- global_sds
  best_subsets <- NULL
  
  # Helper to test a configuration
  test_configuration <- function(current_sds) {
    tryCatch({
      split_dataset_balanced(data, sizes, vars, global_means, current_sds, factors, max_repeats, propThresh = 0.07)
    }, error = function(e) {
      NULL  # Return NULL if the function fails
    })
  }
  
  # Initial attempt with the given global_sds
  best_subsets <- test_configuration(current_sds)
  if (is.null(best_subsets)) {
    stop("No initial solution found with the provided ranges.")
  }
  
  # Iteratively tighten ranges for each variable
  for (var_idx in c(2,1,4,3)) {
    repeat {
      # Tighten the range for the current variable
      current_sds[var_idx] <- current_sds[var_idx] - step
      
      # Test the new configuration
      new_subsets <- test_configuration(current_sds)
      if (!is.null(new_subsets)) {
        # Update the best subsets and continue tightening
        best_subsets <- new_subsets
      } else {
        # Revert to the last successful range and stop tightening this variable
        current_sds[var_idx] <- current_sds[var_idx] + step
        break
      }
    }
  }
  
  return(list(
    subsets = best_subsets,
    final_ranges = current_sds
  ))
}

sceneResult <- optimize_balanced_split(
  data = scenesSelectSet_paired,
  sizes = sizes,
  vars = vars,
  global_means = scene_global_means,
  global_sds = scene_global_sds,
  factors = factors,
  step = 0.05
)


summary(shotSelect3[c(3:8,14)])
summary(shotSelect4_1[c(3:8,14)])
summary(scenesSelect2[c(2:7,13)])

subsetsScenes = sceneResult$subsets

# Check results
for (i in seq(3,1)) {
  cat(paste0("\nSubset ", i, " (n = ", nrow(subsetsScenes[[i]]), "):\n"))
  print("Means:")
  print(colMeans(subsetsScenes[[i]][, vars]))
  print("Factor distributions:")
  for (factor in scene_factors) {
    print(prop.table(table(subsetsScenes[[i]][[factor]])))
  }
}

# Create a named list of subsets
subsets_named <- list("Subset_4s" = subsetsScenes[[3]],
                      "Subset_12s" = subsetsScenes[[2]],
                      "Subset_36s" = subsetsScenes[[1]]
)



# Write to an Excel file
write_xlsx(subsets_named, path = "scene78_subsets02_20241130_02.xlsx")


subsetShots = subsetsScenes
subsetShots[[3]] = shotsSelectSet_paried[1:54,]
subsetShots[[2]] = shotsSelectSet_paried[1:18,]
subsetShots[[1]] = shotsSelectSet_paried[1:6,]

for (i in 3:1) {
  for (ii in 1:length(subsetsScenes[[i]][,1])) {
    whichScene = subsetsScenes[[i]]$whichScene[ii]
    subsetShots[[i]][ii,] = shotsSelectSet_paried[scenesSelectSet_paired$whichScene==whichScene,]
  }
}


for (i in seq(3,1)) {
  cat(paste0("\nSubset ", i, " (n = ", nrow(subsetShots[[i]]), "):\n"))
  print("Means:")
  print(colMeans(subsetShots[[i]][, vars]))
  print("Factor distributions:")
  for (factor in scene_factors) {
    print(prop.table(table(subsetShots[[i]][[factor]])))
  }
}

mean(c(subsetsScenes[[1]]$duration/36,
     subsetsScenes[[2]]$duration/12,
  subsetsScenes[[3]]$duration/4))/
mean(c(subsetShots[[1]]$duration/36,
       subsetShots[[2]]$duration/12,
       subsetShots[[3]]$duration/4))




sceneResult <- optimize_balanced_split(
  data = scenesSelectSet_paired,
  sizes = sizes,
  vars = vars,
  global_means = scene_global_means,
  global_sds = scene_global_sds,
  factors = factors,
  step = 0.05
)




outPath = ".../studyforrest/Scene_subsets_20241205/"
library(stringr)

ratios = c()

for (m in 81:100) {
 outFile = paste0(outPath, "sceneSubset_20241205_", str_pad(m, 3, pad = "0"), ".RData") 
 
 sceneResult_m <- optimize_balanced_split(
   data = scenesSelectSet_paired,
   sizes = sizes,
   vars = vars,
   global_means = scene_global_means,
   global_sds = scene_global_sds,
   factors = factors,
   step = 0.1)
 
 subsetsScenes_m = sceneResult_m$subsets
 
 save(subsetsScenes_m, file = outFile)
 
 subsetShots_m = subsetShots
 for (i in 3:1) {
   for (ii in 1:length(subsetsScenes_m[[i]][,1])) {
     whichScene = subsetsScenes_m[[i]]$whichScene[ii]
     subsetShots_m[[i]][ii,] = shotsSelectSet_paried[scenesSelectSet_paired$whichScene==whichScene,]
   }
 }
 
 ratio = mean(c(subsetsScenes_m[[1]]$duration/36,
                subsetsScenes_m[[2]]$duration/12,
        subsetsScenes_m[[3]]$duration/4))/
   mean(c(subsetShots_m[[1]]$duration/36,
          subsetShots_m[[2]]$duration/12,
          subsetShots_m[[3]]$duration/4))
 
 ratios[m] = ratio
 
}

whichRand = which(ratios == sort(ratios)[2])

inFile = paste0(outPath, "sceneSubset_20241205_", str_pad(whichRand, 3, pad = "0"), ".RData") 
load(inFile)

subsetShots_m = subsetShots

for (i in 3:1) {
  for (ii in 1:length(subsetsScenes_m[[i]][,1])) {
    whichScene = subsetsScenes_m[[i]]$whichScene[ii]
    subsetShots_m[[i]][ii,] = shotsSelectSet_paried[scenesSelectSet_paired$whichScene==whichScene,]
  }
}

# control ratio
mean(c(subsetsScenes_m[[1]]$duration/36,
       subsetsScenes_m[[2]]$duration/12,
       subsetsScenes_m[[3]]$duration/4))/
  mean(c(subsetShots_m[[1]]$duration/36,
         subsetShots_m[[2]]$duration/12,
         subsetShots_m[[3]]$duration/4))



# Check results Scenes
for (i in seq(3,1)) {
  cat(paste0("\nSubset ", i, " (n = ", nrow(subsetsScenes_m[[i]]), "):\n"))
  print("Means:")
  print(colMeans(subsetsScenes_m[[i]][, vars]))
  print("Factor distributions:")
  for (factor in scene_factors) {
    print(prop.table(table(subsetsScenes_m[[i]][[factor]])))
  }
}


# Check results Shots
for (i in seq(3,1)) {
  cat(paste0("\nSubset ", i, " (n = ", nrow(subsetShots_m[[i]]), "):\n"))
  print("Means:")
  print(colMeans(subsetShots_m[[i]][, vars]))
  print("Factor distributions:")
  for (factor in scene_factors) {
    print(prop.table(table(subsetShots_m[[i]][[factor]])))
  }
}


# plot ratio overlap
p = ggplot() +
  geom_histogram(aes(x = c(subsetsScenes_m[[1]]$duration/36,
                           subsetsScenes_m[[2]]$duration/12,
                           subsetsScenes_m[[3]]$duration/4)), color = "darkgreen", fill = "lightgreen", alpha =1, binwidth = 1)+
  geom_histogram(aes(x = c(subsetShots_m[[1]]$duration/36,
                       subsetShots_m[[2]]$duration/12,
                       subsetShots_m[[3]]$duration/4)), color = "darkblue", fill = "lightblue",  alpha =0.5, binwidth = 1) 
p

# Create a named list of subsets
subsets_named <- list("Subset_4s" = subsetsScenes_m[[3]],
                      "Subset_12s" = subsetsScenes_m[[2]],
                      "Subset_36s" = subsetsScenes_m[[1]]
)

# Write to an Excel file
write_xlsx(subsets_named, path = "sceneSubsets_m76_20241205.xlsx")

subsetsShots_named <- list("Subset_4s" = subsetShots_m[[3]],
                      "Subset_12s" = subsetShots_m[[2]],
                      "Subset_36s" = subsetShots_m[[1]]
)

# Write to an Excel file
write_xlsx(subsetsShots_named, path = "shotSubsets_m76_20241205.xlsx")
save.image(".../studyforrest/Scene_classification_newRatios_20241205_final.RData")

p = ggplot() +
  geom_histogram(aes(x = as.numeric(subsets_named[[2]]$whichScene)), color = "darkgreen", fill = "lightgreen", alpha =1, binwidth = 5)+
  geom_histogram(aes(x = as.numeric(subsetsShots_named[[2]]$whichScene)), color = "darkblue", fill = "lightblue",  alpha =0.5, binwidth = 5) 
p



