
function R_spt = VisualComplx_FM(InputImage)

% The Rsuprathreshold metric calculates the visual complexity of an image
% Durmus, D. (2020). Spatial frequency and the performance of image-based visual complexity metrics. IEEE Access.

% based on MATLAB Image Region Analyzes
% https://www.mathworks.com/help/images/calculate-region-properties-using-image-region-analyzer.html

% load an image
if ischar(InputImage) || isstring(InputImage)
    I=imread(InputImage);
elseif isnumeric(InputImage) || islogical(InputImage)
    I = InputImage;
else
    error("Neither Path to image not Image array were provided");
end

[~, ~, numberOfColorChannels] = size(I);

if numberOfColorChannels == 3
    Igray = rgb2gray(I);
else
    Igray = I;
    
end

%convert to binary (graysacale) using 'adaptive' setting
% imbinarize 'global' or 'adaptive'
BinarizedI = imbinarize(Igray,'adaptive');

stats = regionprops(BinarizedI,'Area','Eccentricity');

% % calcualte Rsuprathreshold metric
% total pixel size
RP_tot=size(BinarizedI,1)*size(BinarizedI,2);
% find the number of cells larger than
R_spt=sum([stats.Area]>(RP_tot/25000));

% Meaning of the results
% R_spt < 100 means image is simple
% R_spt > 300 means image is complex
% the rest are medium

end

