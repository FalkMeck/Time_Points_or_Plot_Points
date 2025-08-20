%% MAKE QUEST ORDERS

cd('F:\bigbipsy2\fmecklenbrauck\09_WS24_25\Movie-HINTS_Experiment');

% Within block
% which pictures are flipped? and how many?
block_additionals = cell(size(block_orders)); 
trialPattern = 'trial(\d+)';        % Matches 'trial' followed by one or more digits
scenePattern = '_scene_(\d{3})';   % Matches '_scene_' followed by exactly 3 digits

blocks = {'Scene_4s','Scene_12s', 'Scene_36s','Shot_4s','Shot_12s','Shot_36s'};

for b = 1:numel(blocks)
    content = dir(fullfile(pwd, 'QuestImages', blocks{b})); 
    content = transpose({content(5:end).name}); 
    trials = cellfun(@(x) regexp(x, trialPattern, 'tokens'), content, 'UniformOutput', false);
    trials =cellfun(@(x) x{1}, trials, 'UniformOutput', false);
    trials = cellfun(@(x) str2double(x), trials, 'UniformOutput', false); 
    trials = cell2mat(trials);

    scene = cellfun(@(x) regexp(x, scenePattern, 'tokens'), content, 'UniformOutput', false);
    scene =cellfun(@(x) x{1}, scene, 'UniformOutput', false);
    scene = cellfun(@(x) str2double(x), scene, 'UniformOutput', false); 
    scene = cell2mat(scene);

    quest_blocks.(blocks{b}) = table(content, trials, scene, 'VariableNames', {'pic', 'trial', 'scene'}); 
end

save('quest_blocks_20250109.mat','quest_blocks'); 

%%

sequence = [2, 2, 3, 3, 4, 4];

% Generate all permutations (including duplicates)
allPermutations = perms(sequence);

% Remove duplicate rows to get unique permutations
uniquePermutations = unique(allPermutations, 'rows');


flips = [repmat(2, 16,1); repmat(3, 16,1);repmat(4, 16,1)];

occurance_check = false;
counter = 0; 
while ~occurance_check && counter < 1e6
rand_order= randperm(90,48);
blocked_flips = uniquePermutations(rand_order,:)'; 

    testing = zeros(6,1);
    for i = 1:6
        sorted = sort(blocked_flips(block_orders == i), 'ascend');
        testing(i) = corr(sorted, flips); 
    end
    occurance_check = all(testing > 0.98); 
    counter = counter +1;
end
save('blocked_flips_20250107_03.mat','blocked_flips'); 
for i = 1:48
    disp(i); 
    tabulate(blocked_flips(:,48));
end

for i = 1:6
    disp(i); 
    tabulate(sort(blocked_flips(block_orders == i)));
end


for i = 1:6
    disp(i); 
    tabulate(sort(blocked_flips(i,:)));
end

load('blocked_flips_20250107_03.mat');
load('set_orders_final.mat');
set_orders = set_orders';
load('quest_blocks_20250109.mat');

for b = 1:numel(blocks)
    disp(blocks{b});
    minDeviation = Inf;
    bestVec = [];
    topVec =[];
    for subi = 1:size(set_orders,2)
        how_many_flipped = blocked_flips(b,subi);
        which_fliped = randperm(numel(blocks), how_many_flipped);
        topVec = [topVec, which_fliped];
    end
    % how many flips per condition
    idealFreq = length(topVec) / 6;
    disp(idealFreq);
    
    for iter = 1:1e4
        drawingPool = repmat(1:6,1,ceil(idealFreq));
        vec = [];
        for subi = 1:size(set_orders,2)
            how_many_flipped = blocked_flips(b,subi);
            which_fliped = []; possible = true;
            while length(unique(which_fliped)) ~= how_many_flipped && possible
                which_idx = randperm(length(drawingPool), how_many_flipped);
                which_fliped = drawingPool(which_idx);
                if length(unique(drawingPool)) < how_many_flipped
                    possible = false;
                end
            end
            
            vec = [vec, which_fliped];
            for num = which_fliped
                % Find the index of the first occurrence of the drawn number
                idx = find(drawingPool == num, 1);
                % Remove the first occurrence of that number
                drawingPool(idx) = [];
            end
        end
        
        if length(vec) == length(topVec)
            counts = histcounts(vec , 0.5:6.5);
            % Calculate the deviation from the ideal frequency
            
            deviation = sum(abs(counts - idealFreq));
            
            % Update the best permutation if this one is better
            if deviation < minDeviation
                minDeviation = deviation;
                bestVec = vec;
            end
            
            % Break early if a perfect match is found
            if deviation == 0 && possible
                disp('found perfection');
                break;
            end
        end
    end
    disp(counts); disp(deviation); disp(possible); 
    flips.(blocks{b}) = bestVec;
end

%% Export Flips
for b = 1:numel(blocks)
    disp(blocks{b}); 
    blockDir = fullfile(pwd, 'quest_orders', blocks{b});
    if ~exist(blockDir, 'dir')
        mkdir(blockDir);
    end
    which_fliped_all = correct_flips4.(blocks{b});
    startP = 0; endP = 0;
    testing = zeros(1,6); testing_in_order = zeros(1,6); 
     counts = histcounts(which_fliped_all , 0.5:6.5);
     disp(counts); 
    for subi = 1:size(set_orders,2)
        how_many_flipped = blocked_flips(b,subi);
        endP = startP+how_many_flipped;
        startP = startP+1; 
        %disp(["how many flips: ", num2str(how_many_flipped),...
         %   ", Range: ", num2str(endP-startP+1)]); 
        which_fliped = which_fliped_all(startP:endP);
        if how_many_flipped ~= length(unique(which_fliped))
             disp(subi);
             disp(which_fliped);
         end        
        outTable = quest_blocks.(blocks{b});
        outTable.seen = ones(6,1);
        outTable.seen(which_fliped) = repmat(2, how_many_flipped, 1);
        testing(which_fliped) = testing(which_fliped) +1;
        outTableSorted = sortrows(outTable,2,'ascend');
        testing_in_order(outTableSorted.seen == 2) = testing_in_order(outTableSorted.seen == 2) +1;
        outName = [blockDir,filesep,'quest_input_sub',sprintf('%02d',subi),'.txt'];
       % writetable(outTableSorted, outName, 'delimiter', '\t');
        startP = endP;
    end
    disp(testing);
    disp(testing_in_order); 
end
