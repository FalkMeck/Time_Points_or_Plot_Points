#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:44:17 2024

@author: f_meck01
"""
# IMPORT LIBRARIES
import pandas as pd
import subprocess


# DEFINE DIRECTORIES
movieDir = .../Movie/"
outDir = movieDir + "_Scenes_Shots_v1/"



# Load CSV with start and end times in HH:MM:SS:FF format
df = pd.read_csv(movieDir + "all4extract_???.csv") # Filre that contains both selected scenes and shots
df['place_cleaned'] = df['place'].str.replace(r'[^a-zA-Z]', '_', regex=True)

# Specify input video file path
input_video = movieDir + "ForrestGump_StudyCut.mp4" # Specific movie Cut (ATTENTION: StudyForrest used a version that excluded all scene with discriminating slurs) 

# Function to convert HH:MM:SS:FF to HH:MM:SS.sss
def convert_to_hhmmss_ms(timecode, fps=25):
    hh, mm, ss, ff = map(int, timecode.split(':'))
    # Convert frame count to milliseconds
    milliseconds = int((ff / fps) * 1000)
    return f"{hh:02}:{mm:02}:{ss:02}.{milliseconds:03}"

# Iterate over each scene and cut the segment
for index in range(len(df)):
    start = convert_to_hhmmss_ms(df['start_time_buffer'][index])  # Convert start time
    end = convert_to_hhmmss_ms(df['end_time_buffer'][index])      # Convert end time
    
    file_name = f"{df['hierarchy'][index]}_{df['duration'][index]:02}_scene_{df['whichScene'][index]:03}_{df['place_cleaned'][index]}_{df['time'][index]}_{df['location'][index]}.mp4"
    output_file = outDir + file_name

    # Run FFmpeg command to cut the segment with 25 FPS interpretation
    command = [
        "ffmpeg", "-i", input_video, "-ss", start, "-to", end,
        "-c", "copy", output_file
    ]
    
    subprocess.run(command)
    print(f"Created {output_file} from {start} to {end}.")