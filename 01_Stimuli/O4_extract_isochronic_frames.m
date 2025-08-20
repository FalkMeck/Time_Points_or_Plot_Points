%% EXTRACT ISOCHRONOUS FRAMES FROM VIDEO CLIPS
cd('...\Movie'); 

clipInfo = readtable('all4extract_m76_20241205.csv');
clipInfo.place_cleaned = regexprep(clipInfo.place, '[^a-zA-Z]', '_');

fps = 3;
clipInfo.duration_buffer_ratio = clipInfo.vid_duration_buffer./(clipInfo.duration * fps);

for i = 1:size(clipInfo,1)
    filename = [clipInfo.hierarchy{i}, '_', ...
            sprintf('%02d', clipInfo.duration(i)), '_scene_', ...
            sprintf('%03d', clipInfo.whichScene(i)), '_',...
            clipInfo.place_cleaned{i},'_',...
            clipInfo.time{i},'_',clipInfo.location{i}];
    filePath = fullfile(pwd, '_Scenes_Shots_v1', [filename, '.mp4']); 
    
    numberExtract = clipInfo.duration(i)*fps;
    
    extract_iso_frames(filePath, numberExtract);
    disp(['DONE: ', num2str(i), '/', num2str(size(clipInfo,1))]); 
end

function [] = extract_iso_frames(vidPath, frames2extract)
    
    [vidFolder,name,~] = fileparts(vidPath); 
    disp(['Frame extraction for: ', name]);
    % read video object
    vidObj = VideoReader(vidPath);
    vidframes = read(vidObj,[1 Inf]);
    clear vidObj;
    
    % shouldn't happen, than 
    if frames2extract >= size(vidframes,4)
        disp([vidPath, ': length did not work']);
        number_of_frames = frames2extract;
    else
        number_of_frames = size(vidframes,4);
    end
    
    steps = number_of_frames/frames2extract;
    whichFrames = round(1:steps:number_of_frames,0);
    % just control, if it is perfectly dividable, there might be an issue
    if (numel(whichFrames)-1) == frames2extract
        whichFrames = whichFrames(2:end);
    elseif numel(whichFrames) == frames2extract
        whichFrames = whichFrames;
    else
       disp([vidPath, ': whichframes did not work']);      
    end
    
    % set outDirectory
    outDir = fullfile(vidFolder, name);
    if ~exist(outDir, 'dir')
        mkdir(outDir); 
    end
    
    % for all frames, save images 
    fCount = 0; 
    for ii = whichFrames
        fCount = fCount +1; 
        disp(['exporting frame ',sprintf('%03d', fCount),': ', sprintf('%04d', ii)]); 
        frame = vidframes(:,:,:,ii); 
        frameName = fullfile(outDir, [name, '_', sprintf('%04d', ii),'.png']);
        imwrite(frame, frameName);
    end

end