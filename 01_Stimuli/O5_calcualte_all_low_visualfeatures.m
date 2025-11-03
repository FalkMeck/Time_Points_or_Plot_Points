%% calcualte low-level visual features: luminace, contrast, color saturation, (motion energy) and overall visual complexity

cd('...');
out_dir = fullfile(pwd, '_Visual_Features');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
addpath('...\01_Stimuli\AdditionalFiles');

imgRoot = '...\FrameStimuli_cropped';
ContLvl = {'Scene','Shot'};
Dur = {'4s', '12s', '36s'};

%%
cols = {'ContLvl', 'Duration', 'trial', 'place','time','location','img',...
    'lum', 'contrast', 'saturation', 'complexity'};


outCell = cell(648*6,numel(cols));

m = 0;
for h = 1:numel(ContLvl)
    disp(ContLvl{h});
    for d = 1:numel(Dur)
        disp(Dur{d});
        
        imgDir = fullfile(imgRoot, [ContLvl{h},'_',Dur{d}]);
        imgs = dir(fullfile(imgDir,'*.png*'));
        imgs = {imgs.name};
        
        for i = 1:numel(imgs)
            m = m+1;
            % read image
            outCell(m,1) = ContLvl(h); outCell(m,2) = Dur(d);
            imgSplit = strsplit(imgs{i},'_');
            
            outCell{m,3} = str2double(extract(imgSplit{1}, digitsPattern));
            outCell{m,7} = str2double(extract(imgSplit{end}, digitsPattern));
            outCell(m,6) = imgSplit(end-1);
            outCell(m,5) = imgSplit(end-2);
            outCell{m,4} = strjoin(imgSplit(6:(end-3)),'_');
            imageArray = imread(fullfile(imgDir,imgs{i}));
            
            [~, ~, numberOfColorChannels] = size(imageArray);
            
            %% luminance
            % https://de.mathworks.com/matlabcentral/answers/858825-how-to-get-luminance-of-an-image
            % based on gray scale image intensities
            %% contrast
            %https://de.mathworks.com/matlabcentral/answers/595216-how-to-get-an-image-contrast-value
            % difference between minimal and maximal gray level intensity
            if numberOfColorChannels == 3
                % luminace
                outCell{m,8} = mean2(rgb2gray(imageArray));
                % contrast
                grayImage = double(rgb2gray(imageArray));
                %outCell{m,9} = max(grayImage(:)) - min(grayImage(:));
                %RMS contrast measures intensity variability across pixels, more robust for natural images:
                outCell{m,9} = std(grayImage(:)) / mean(grayImage(:));
            else
                % luminace
                outCell{m,8} = mean2(imageArray);
                % contrast
                grayImage = double(imageArray);
                %outCell{m,9} = max(grayImage(:)) - min(grayImage(:));
                %RMS contrast measures intensity variability across pixels, more robust for natural images:
                outCell{m,9} = std(grayImage(:)) / mean(grayImage(:));
            end
            %% Saturation
            %https://de.mathworks.com/matlabcentral/answers/65000-how-can-i-calculate-brightness-contrast-hue-and-saturation-of-a-image
            imageArray_HSV = rgb2hsv(imageArray);
            outCell{m,10} = mean2(imageArray_HSV(:,:,2));
            
            %% overall complexity
            %https://de.mathworks.com/matlabcentral/fileexchange/107924-visual-complexity-metric-for-images
            % Durmus, D. (2020). Spatial Frequency and the Performance of Image-Based Visual Complexity Metrics. IEEE Access, 8, 100111–100119. Institute of Electrical and Electronics Engineers (IEEE). Retrieved from https://doi.org/10.1109%2Faccess.2020.2998292
            outCell{m,11} = VisualComplx_FM(imageArray);
            
            
            disp(['Done: ', num2str(m), '/', num2str(size(outCell,1))]);
            
        end
    end
end

T = cell2table(outCell);
T.Properties.VariableNames = cols;
writetable(T, fullfile(out_dir, 'LowVisualFeature.csv'));


%% motion energy per video frame sequence

cols = {'ContLvl', 'Duration', 'trial', 'place','time','location','img','motion_energy'};
outCell = cell(156,numel(cols));
m = 0;

for h = 1:numel(ContLvl)
    disp(ContLvl{h});
    for d = 1:numel(Dur)
        disp(Dur{d});
        
        imgDir = fullfile(imgRoot, [ContLvl{h},'_',Dur{d}]);
        imgs = dir(fullfile(imgDir,'*.png*'));
        imgs = {imgs.name};
        
        imgSplits = cellfun(@(x) strsplit(x,'_'), imgs, 'UniformOutput', false);
        imgTrials = cellfun(@(x)str2double(extract(x{1}, digitsPattern)), imgSplits, 'UniformOutput', false);
        
        imgTrials = cell2mat(imgTrials);
        blockTrials = unique(imgTrials);
        
        for b = 1:length(blockTrials)
            m = m+1;
            imageSubList = imgs(imgTrials == blockTrials(b));
            
            outCell(m,1) = ContLvl(h); outCell(m,2) = Dur(d);
            imgSplit = strsplit(imageSubList{1},'_');
            
            outCell{m,3} = str2double(extract(imgSplit{1}, digitsPattern));
            outCell{m,7} = str2double(extract(imgSplit{end}, digitsPattern));
            outCell(m,6) = imgSplit(end-1);
            outCell(m,5) = imgSplit(end-2);
            outCell{m,4} = strjoin(imgSplit(6:(end-3)),'_');
            
            imageArray = imread(fullfile(imgDir,imageSubList{1}));
            
            I = zeros(size(imageArray,1), size(imageArray,2), numel(imageSubList)); 
            
            for i = 1:numel(imageSubList)
                imgTmp = imread(fullfile(imgDir,imageSubList{i}));
                I(:,:,i) = im2double(rgb2gray(imgTmp)); 
            end
            
            motion_energy = zeros(1, size(I,3)-1);
            
            for t = 2:size(I,3)
                frame_prev = I(:,:,t-1);
                frame_curr = I(:,:,t);
                motion_energy(t-1) =  mean((frame_curr(:) - frame_prev(:)).^2);% or sum
            end          
            
            %average motion energy per trial
            outCell{m,8} = mean(motion_energy); 
            disp(['Done: ', num2str(m), '/', num2str(size(outCell,1))]);
            
        end
        
        
    end
end

T = cell2table(outCell);
T.Properties.VariableNames = cols;
writetable(T, fullfile(out_dir, 'Motion_energy.csv'));