# Header
scenario = "Movie_HINTS.sce";
pcl_file = "Movie_HINTS.pcl";

default_font_size = 8;
default_font = "Microsoft JhengHei Light";	
default_background_color = 20,20,20;###oder:150,150,150
default_text_color = 230, 230, 230;###oder:255,255,255

screen_width=1920;
screen_height=1080;
screen_bit_depth = 32;

response_matching = simple_matching;	################
response_logging = log_all;				################

active_buttons = 3;      
button_codes = 1,2,3;

#Fontgroessen
$font_fixcross = 135;
$font_delimiter = 68;
$font_text = 40;
$font_numbers = 180;

# Koordinaten
$x_coord_zentral = 0;
$y_coord_zentral = 0; 

# Header End

#scenario_type = fMRI;  # für Scanner
scenario_type = fMRI_emulation;  # für Scanner
scan_period = 1500;						################
pulses_per_scan = 1; 					################
pulse_code = 55; 							################


begin;

##########################################################################
#		
#		Text
#
##########################################################################

## Zeichen
text {caption = "+"; system_memory = true; font_size = $font_fixcross; } fixcross_text;
text {caption = " "; system_memory = true; font_size = $font_delimiter; } delimiter_text;

text {caption = "Willkommen!

Dieses Experiment besteht aus zwei Teilen: 
einer Resting-State-Messung und dem Hauptteil, 
bei dem Du eine Aufgabe bekommst.

Wir beginnen mit der Resting-State-Messung:
Bleibe dazu bitte einfach entspannt liegen, 
bewege Dich nicht und 
betrachte das Fixationskreuz

Wenn Du bereit bist, 
drücke bitte auf eine der beiden Antwort-Tasten"; font_size = $font_text;} Rest_Start_text;

text {caption = "Einen kleinen Moment noch,
es geht sofort weiter..."; font_size = $font_text;} Rest_Wait_text;

text {caption = "Die Resting-State-Messung ist geschafft! 

Im Hauptteil des Experiments siehst Du gleich 
6 Blöcke von Frames aus dem Film 'Forrest Gump'. 
Betrachte die Frames aufmerksam, als würdest Du den Film schauen.

Nach jedem Block siehst Du 6 Frames aus dem vorherigen Block. 
Die Frames können WIE IM BLOCK oder GESPIEGELT präsentiert.
Bewerte, ob die Frames WIE IM BLOCK gezeigt werden. 

ACHTUNG: Die Zuordnung der Antwort-Tasten kann wechseln.

Alles klar? Wenn Du bereit bist,
drücke bitte auf eine der beiden Antwort-Tasten."; font_size = $font_text;} Start_text; # change that it starts automatically?

text {caption = "Jetzt geht es mit der Aufgabe weiter.

Werden die Frames WIE IM BLOCK präsentiert?"; font_size = $font_text;} Before_task_text;

text {caption = "Super!\n Du hast es geschafft!

Bleibe bitte trotzdem noch ruhig liegen."; font_size = $font_text;} End_text;

text {caption = " "; font_size = $font_text;} break_text;

text {caption = "Kurze Pause...

Gleich geht es mit"; font_size = $font_text;} next_block_text_break;
text {caption = "Bereit?

Gleich startet"; font_size = $font_text;} next_block_text_start;
text {caption = "Block N"; font_size = $font_text;} next_block_text_blk;
text {caption = "weiter."; font_size = $font_text;} next_block_text_weiter;

text {caption = "Korrekt eingeschätzte Frames:"; font_size = $font_text;} final_feedback_text;
text {caption = "0/36"; font_size = $font_text;} final_score_text;



text {caption = "WIE IM BLOCK?"; font_size = $font_text;} question_text;
text {caption = "links"; font_size = $font_text;} left_text;
text {caption = "rechts"; font_size = $font_text;} right_text; 

#text {caption = "Richtig"; font_size = $font_text; font_color = 0, 255, 0;} feedback_true_text;
#text {caption = "Falsch"; font_size = $font_text; font_color = 255, 0, 0;} feedback_false_text;
#text {caption = "Nicht gedrückt"; font_size = $font_text; font_color = 230, 230, 230;} feedback_miss_text;

#text {caption = "Loading frames..."; font_size = $font_text;} loading_text;

##########################################################################
#		
#		Frames
#
##########################################################################

bitmap { filename = "Marius.png"; scale_factor = 1.5;} frame;
bitmap { filename = "Marius.png"; scale_factor = 1.5;} quest_frame;

##########################################################################
#		
#		Pictures
#
##########################################################################

picture {} default;

# define instruction, fixation, picture and video
picture { text Rest_Start_text; x = 0; y = 0;} Rest_Start_pic;

picture { text Rest_Wait_text; x = 0; y = 0;} Rest_Wait_pic;

picture { text Start_text; x = 0; y = 0; } Start_pic;

picture { text fixcross_text; x = $x_coord_zentral; y = $y_coord_zentral; } fixcross_pic;

picture { text Before_task_text; x=0; y=0;} Before_task_pic;

picture { text break_text; x = 0; y = 0;} break_pic;

picture { text End_text; x=0; y=0;} End_pic;

picture {bitmap frame; x=0; y=0;} frame_pic;

picture {bitmap quest_frame; x=0; y=0;} quest_frame_pic;

picture{
	text question_text; x = 0; y = 0;
	text left_text; x = -600; y = 0;
	text right_text; x = 600; y = 0;} question_pic;
	
# picture{text feedback_true_text; x = 0; y = 0;} feedback_pic;

picture{text final_feedback_text; x = 0; y = 300;
		  text final_score_text; x = 0; y = 0;} final_feedback_pic;
		
picture{text next_block_text_break; x = 0; y = 150;
		  text next_block_text_blk; x = 0; y = 0;
		  text next_block_text_weiter; x = 0; y = -150;} block_break_pic;

picture{text next_block_text_start; x = 0; y = 75;
		  text next_block_text_blk; x = 0; y = -75;} block_start_pic;

#picture{text loading_text; x = 0; y = 0;} loading_pic;

##########################################################################
#		
#		Trials (also includes some unused ones)
#
##########################################################################

# Welcome and Ready trial
trial {
	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 1,2;
	picture Rest_Start_pic;
} Rest_Start_trial;


trial {
	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 3;
	picture Rest_Wait_pic;
} Rest_Wait_trial;


trial {
	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 1,2;
	picture Start_pic;
} Start_trial;


#Frame Trial
trial {
   trial_type = fixed;
	trial_duration = 333;
	clear_active_stimuli = false;
	   stimulus_event {
		picture frame_pic;
		time = 0;
		code = "Frame";
		response_active = false;
   } frame_event;
  }frame_trial; 

# Before task trial
trial{trial_type = fixed;
		trial_duration = 1500;
		picture Before_task_pic;}Before_task_trial;

#Quest Frame Trial
trial {
   trial_type = fixed;
	trial_duration = 4500;
	clear_active_stimuli = false;
	   stimulus_event {
		picture quest_frame_pic;
		time = 0;
		code = "Quest_Frame";
		response_active = false;
   } quest_frame_event;
  }quest_frame_trial; 

# Question Trial
trial {
	trial_type = specific_response;
	terminator_button = 1,2; 
	trial_duration = 1500; 
	picture default; 
	stimulus_event {
		picture question_pic;
		time = 0;
		code = "Quest";
	} question_event;
} question_trial;

# Feedback Trial
#trial{
#	trial_type = fixed; 
#	trial_duration = 500;
#	picture feedback_pic;} feedback_trial; 

# Fixcross Trial
trial {
   trial_type = fixed;
	trial_duration = 500;
		stimulus_event {	
			picture fixcross_pic;
			code="Fix";
			time=0;  
		} fix_event;
} fixcross_trial;

# Fixcorss blink trial 
trial{
	trial_type = fixed;
	trial_duration = 333;
	stimulus_event{
		picture fixcross_pic;
		code ="FixBlink";
		time = 0;
		duration = 166;}fixblink_event;} fixblink_trial;

# Break
trial{trial_type = fixed;
		trial_duration = 100;
		picture break_pic;}break_trial;
		
# block break
trial{trial_type = fixed; 
		trial_duration = 24;
		picture block_break_pic;}block_break_trial;

# block start coordination		
trial{trial_type = specific_response; 
		terminator_button = 3; #wird extern vom Kontrollraum weiter gedrückt mit E 
		trial_duration = forever;
		picture block_start_pic;}block_start_trial;

# loading wait trial
trial{trial_type = fixed; 
		trial_duration = 24;
		picture Rest_Wait_pic;}loading_trial;
		
		
# Final Feedback Trial
trial{trial_type = fixed;
		trial_duration = 2500;
		picture final_feedback_pic;}final_feedback_trial;
		
#EndeTrial
trial { 
			trial_type= specific_response;
			trial_duration = forever;
			terminator_button = 3;
			picture End_pic;} End_trial;
			
			
# progress bar
# borders
#box { height = 2; width = 1500; color = 230,230,230; } horiz;
#box { height = 100; width = 2; color = 230,230,230; } vert;
#box { height = 100; width = 1500; color = 230,230,230;} bar; 


#trial {
	#trial_duration = 100;
	#picture{
	#	text{caption = "caption"; font_size = $font_text;}short_break; x = 0; y = 0;
   #   box horiz; x = 0; y = -350;
   #   box horiz; x = 0; y = -450;
   #   box vert; x = -750; y = -400;
   #   box vert; x = 750; y = -400;
#		box bar; x = 0; y = -400;}bar_pic;
#}progress_bar;

#trial {trial_type = specific_response; terminator_button = 3; #wird extern vom Kontrollraum weiter gedrückt mit E#
#		trial_duration = forever;
#		picture{text{caption = "caption"; font_size = $font_text;}short_break_cont; x = 0; y = 0;}short_break_pic;
#		code = "break";
#} br_trial;



###################################################
###				   fMRT Standard			   	   ### #????		################
###################################################

# for write_log subroutine	#das ist ein Standard trial - IMMER SO UEBERNEHMEN BEI MR-MESSUNGEN
#trial {	
#stimulus_event {
#nothing {};
#time = 0;
#code = "nothing";	#diesen code machen wir uns nacher zu nutze
#}evt_log;
#} param_log;

