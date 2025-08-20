# Balance the step order

setwd(".../studyforrest/gump_emotions")
library(stringr)
library(ggplot2)
## LOAD SCENES
# Get the names of all sheets
sheet_names <- readxl::excel_sheets("sceneSubsets_???.xlsx")

# Read each sheet into a list of data frames
scenes <- lapply(sheet_names, function(sheet) {
  readxl::read_excel("sceneSubsets_???.xlsx", sheet = sheet)
})

## LOAD SHOTS
sheet_names <- readxl::excel_sheets("shotSubsets_???.xlsx")

# Read each sheet into a list of data frames
shots <- lapply(sheet_names, function(sheet) {
  readxl::read_excel("shotSubsets_???.xlsx", sheet = sheet)
})


shots4select = readxl::read_excel("./AddtionalFiles/shots4select_corrected.xlsx")
scenes4select = readxl::read_excel("./AddtionalFiles/scenes4select.xlsx")
cumSumShots = readxl::read_excel("./AddtionalFiles/cumSumShots.xlsx")
names(cumSumShots) = "cumSumShot"
cumSumShots[68:205,1] = cumSumShots[68:205,1] -1

scenes4select$start_time_shots = shots4select$start_time[cumSumShots$cumSumShot[1:205]]
scenes4select$start_time_shots[68] = scenes4select$start_time[68] # correct WW2
shot_timecodes_last = "01:58:02:05"
scenes4select$end_time_shots = c(scenes4select$start_time_shots[2:length(scenes4select$start_time)], shot_timecodes_last)

# Correct for WWII
shots[[1]]$whichScene[shots[[1]]$whichScene == 68] = 69

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
}, scenes4select$start_time, scenes4select$start_time_shots)

# Print the differences
# frame_differences

shots4select$end_time = c(shots4select$start_time[2:length(shots4select$start_time)], shot_timecodes_last)


## create data frame with all shots/scenes to cut
# Function to convert HH:MM:SS:FF to total frames
convert_to_frames <- function(time, fps) {
  # Split the time into components
  time_parts <- strsplit(time, ":")[[1]]
  hours <- as.numeric(time_parts[1])
  minutes <- as.numeric(time_parts[2])
  seconds <- as.numeric(time_parts[3])
  frames <- as.numeric(time_parts[4])
  
  # Convert to total frames
  total_frames <- frames + (seconds * fps) + (minutes * 60 * fps) + (hours * 3600 * fps)
  return(total_frames)
}

scenes4extract = data.frame(whichScene = rep(0,78),
                            hierarchy = rep("Scene",78),
                            duration = 0,
                            place = NA,
                            time = NA,
                            location = NA,
                            start_time = NA,
                            end_time = NA,
                            vid_duration = NA,
                            howManyShots = 0,
                            arousal_mean = 0,
                            number_people = 0)

m = 1
durations = c(4,12,36)
for (i in 1:3) {
  for (ii in 1:length(scenes[[i]]$whichScene)) {
    scenes4extract$whichScene[m] = as.numeric(scenes[[i]]$whichScene[ii]) # as.numeric(row.names(scenes4select)[scenes4select$start_time_shots == scenes[[i]]$start_time_shots[ii]]) # which Scene
    scenes4extract$duration[m] = durations[i] # duration
    scenes4extract$place[m] = scenes[[i]]$place[ii] # which Place
    scenes4extract$time[m] = scenes[[i]]$time[ii] # which Time
    scenes4extract$location[m] = scenes[[i]]$location[ii] # which Place
    scenes4extract$start_time[m] = scenes4select$start_time_shots[scenes4extract$whichScene[m]] # which start_time according to the shots
    scenes4extract$end_time[m] = scenes4select$end_time_shots[scenes4extract$whichScene[m]] # which end_time according to the shots
    scenes4extract$vid_duration[m] = convert_to_frames(scenes4extract$end_time[m], 25) - 
      convert_to_frames(scenes4extract$start_time[m], 25) # calculate duration in frames
    scenes4extract$howManyShots[m] = scenes[[i]]$howManyShots[ii] # how many shots
    scenes4extract$arousal_mean[m] = scenes[[i]]$arousal_mean[ii] # arousal
    scenes4extract$number_people[m] = scenes[[i]]$peoplePerScene[ii] # people per scene
    m = m+1
  }
}



# and for shots
shots4extract = data.frame(whichScene = rep(0,78),
                           hierarchy = rep("Shot",78),
                           duration = 0,
                           place = NA,
                           time = NA,
                           location = NA,
                           start_time = NA,
                           end_time = NA,
                           vid_duration = NA,
                           howManyShots = 1,
                           arousal_mean = 0,
                           number_people = 0)

m = 1
durations = c(4,12,36)
for (i in 1:3) {
  for (ii in 1:length(shots[[i]]$start_time)) {
    shots4extract$whichScene[m] = shots[[i]]$whichScene[ii] # which Scene
    shots4extract$duration[m] = durations[i] # duration
    shots4extract$place[m] = scenes4select$place[shots4extract$whichScene[m]] # which Place
    shots4extract$time[m] = scenes4select$time[shots4extract$whichScene[m]] # which Time
    shots4extract$location[m] = scenes4select$location[shots4extract$whichScene[m]] # which Place
    shots4extract$start_time[m] = shots4select$start_time[shots4select$start_time == shots[[i]]$start_time[ii]] # which start_time according to the shots
    shots4extract$end_time[m] = shots4select$end_time[shots4select$start_time == shots[[i]]$start_time[ii]] # which end_time according to the shots
    shots4extract$vid_duration[m] = convert_to_frames(shots4extract$end_time[m], 25) - 
      convert_to_frames(shots4extract$start_time[m], 25) # calculate number of frames
    shots4extract$arousal_mean[m] = shots[[i]]$arousal_mean[ii] # arousal
    shots4extract$number_people[m] = shots[[i]]$peoplePerScene[ii]# people per scene
    m = m+1
  }
}

scenes4extract$duration_ratio = scenes4extract$vid_duration/(scenes4extract$duration * 3)
shots4extract$duration_ratio = shots4extract$vid_duration/(shots4extract$duration * 3)

tapply(scenes4extract$duration_ratio, scenes4extract$duration, mean)/
  tapply(shots4extract$duration_ratio, shots4extract$duration, mean)


p = ggplot() +
  geom_histogram(data = scenes4extract, aes(x = vid_duration), color = "darkgreen", fill = "lightgreen", alpha =1, binwidth = 100)+
  geom_histogram(data = shots4extract, aes(x = vid_duration), color = "darkblue", fill = "lightblue",  alpha =0.5, binwidth = 100) 
p

p = ggplot() +
  geom_histogram(data = scenes4extract, aes(x = arousal_mean), color = "darkgreen", fill = "lightgreen", alpha =1, binwidth = 0.1)+
  geom_histogram(data = shots4extract, aes(x = arousal_mean), color = "darkblue", fill = "lightblue",  alpha =0.5, binwidth = 0.1) 
p

p = ggplot() +
  geom_histogram(data = scenes4extract, aes(x = duration_ratio), color = "darkgreen", fill = "lightgreen", alpha =1, binwidth = 10)+
  geom_histogram(data = shots4extract, aes(x = duration_ratio), color = "darkblue", fill = "lightblue",  alpha =0.5, binwidth = 10) 
p


# combine both data sets

all4extract = rbind(scenes4extract, shots4extract)

# create buffer if cuts are weird
parts = str_split(all4extract$start_time, pattern = ":")

all4extract$start_time_buffer = NA

for (i in 1:length(parts)) {
  hours = as.integer(parts[[i]][1])
  minutes = as.integer(parts[[i]][2])
  seconds = as.integer(parts[[i]][3])
  frames = as.integer(parts[[i]][4])
  
  # Add 4 frames
  frames <- frames + 4
  
  # Adjust frames and seconds if frames exceed 24
  if (frames >= 25) {
    frames <- frames - 25
    seconds <- seconds + 1
  }
  
  # Adjust seconds and minutes if seconds exceed 59
  if (seconds >= 60) {
    seconds <- seconds - 60
    minutes <- minutes + 1
  }
  
  # Adjust minutes and hours if minutes exceed 59
  if (minutes >= 60) {
    minutes <- minutes - 60
    hours <- hours + 1
  }
  
  # Reformat the new timecode as HH:MM:SS:FF
  all4extract$start_time_buffer[i] = sprintf("%02d:%02d:%02d:%02d", hours, minutes, seconds, frames)
}

parts = str_split(all4extract$end_time, pattern = ":") 
all4extract$end_time_buffer = NA

for (i in 1:length(parts)) {
  hours = as.integer(parts[[i]][1])
  minutes = as.integer(parts[[i]][2])
  seconds = as.integer(parts[[i]][3])
  frames = as.integer(parts[[i]][4])
  
  # Add 4 frames
  frames <- frames - 4
  
  # Adjust frames and seconds if frames exceed 24
  if (frames < 0 ) {
    frames <- frames + 25
    seconds <- seconds - 1
  }
  
  # Adjust seconds and minutes if seconds exceed 59
  if (seconds < 0) {
    seconds <- seconds + 60
    minutes <- minutes - 1
  }
  
  # Adjust minutes and hours if minutes exceed 59
  if (minutes < 0) {
    minutes <- minutes + 60
    hours <- hours - 1
  }
  
  if (hours < 0) {
    hours <- 0
    minutes <- 0
    seconds <- 0
    frames <- 0
  }
  
  # Reformat the new timecode as HH:MM:SS:FF
  all4extract$end_time_buffer[i] = sprintf("%02d:%02d:%02d:%02d", hours, minutes, seconds, frames)
}

all4extract$vid_duration_buffer = all4extract$vid_duration -8


# BALANCE ORDER ####

list4extract = list(all4extract[all4extract$hierarchy == "Scene" & all4extract$duration == 4,], 
                    all4extract[all4extract$hierarchy == "Scene" & all4extract$duration == 12,], 
                     all4extract[all4extract$hierarchy == "Scene" & all4extract$duration == 36,], 
                    all4extract[all4extract$hierarchy == "Shot" & all4extract$duration == 4,], 
                    all4extract[all4extract$hierarchy == "Shot" & all4extract$duration == 12,], 
                    all4extract[all4extract$hierarchy == "Shot" & all4extract$duration == 36,])

balance_clip_order <- function(blocks, timestamp_col, location_col, scene_col, fps = 25, max_repeats = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Helper function to convert timestamps to total frames
  timestamp_to_frames <- function(timestamp, fps) {
    sapply(timestamp, function(t) {
      parts <- as.numeric(strsplit(t, ":")[[1]])
      frames <- parts[1] * 3600 * fps + parts[2] * 60 * fps + parts[3] * fps + parts[4]
      return(frames)
    })
  }
  
  # Balance a single block
  balance_single_block <- function(block, timestamp_col, location_col, scene_col, fps, max_repeats) {
    timestamps <- block[[timestamp_col]]
    block[["frame_position"]] <- timestamp_to_frames(timestamps, fps)
    
    best_order <- NULL
    best_score <- Inf
    
    for (i in 1:max_repeats) {
      order <- sample(1:nrow(block))  # Randomly permute order
      reordered_block <- block[order, ]
      
      # Calculate jump distances
      reordered_frames <- block[["frame_position"]][order]
      jumps <- diff(reordered_frames)
      forward_jumps <- jumps[jumps > 0]
      backward_jumps <- jumps[jumps < 0]
      jump_balance_score <- var(c(sum(forward_jumps), abs(sum(backward_jumps))))
      # print(jump_balance_score)
      # Check location constraint
      location_valid <- all(reordered_block[[location_col]][-1] != reordered_block[[location_col]][-nrow(reordered_block)])
      
      # Check scene constraint
      scene_valid <- all(reordered_block[[scene_col]][-1] != reordered_block[[scene_col]][-nrow(reordered_block)])
      
      # Compute score only if constraints are satisfied
      if (location_valid && scene_valid) {
        score <- jump_balance_score
        if (score < best_score) {
          best_order <- order
          best_score <- score
        }
      }
    }
    
    if (is.null(best_order)) {
      stop("Could not find a valid order within max_repeats.")
    }
    
    block = block[best_order, ]
    block[["Jumps"]] = c(NA,diff(block[["frame_position"]]))
    return(block)
  }
  
  # Apply balancing to each block
  balanced_blocks <- lapply(blocks, function(block) {
    balance_single_block(block, timestamp_col, location_col, scene_col, fps, max_repeats)
  })
  
  return(balanced_blocks)
}

balanced_blocks_list = balance_clip_order(list4extract,
                                          timestamp_col = "start_time_buffer",
                                          location_col = "place",
                                          scene_col = "whichScene",
                                          fps = 25,
                                          max_repeats = 1e5,
                                          seed = 420)
View(balanced_blocks_list[[1]])

save(balanced_blocks_list, file = "balanced_blocks_list.RData")
balancedBlocks4export = rbind(balanced_blocks_list[[1]],
                              balanced_blocks_list[[2]],
                              balanced_blocks_list[[3]],
                              balanced_blocks_list[[4]],
                              balanced_blocks_list[[5]],
                              balanced_blocks_list[[6]])
write.csv(balancedBlocks4export, file = "balancedBlocks4export.csv", row.names = F)
