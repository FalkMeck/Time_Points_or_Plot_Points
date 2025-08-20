%% Block order Movie-HINTS

% the experiment is going to have 6 Blocks of teh following conditions:
% 1. Scene 4 seconds
% 2. Scene 12 seconds
% 3. Scene 36 seconds
% 4. Shot 4 seconds
% 5. Shot 12 seconds
% 6. shot 36 seconds

% Across participants we want to balance the start condtion, the
% transistions between condtions, possibly also how often each condition is
% in each position in the experiment once

%% genereate iteratively

cd('...\Movie');

nBlocks = 6;
max_retries = 100;
max_rounds = 8;

perfect_matrix = ones(nBlocks, nBlocks);
perfect_matrix = perfect_matrix - diag(diag(perfect_matrix));

max_for = 1e5;

for file = 1:20
    isDone = false;
for m = 1:max_for
    round = 0;
    disp([num2str(m), '/', num2str(max_for)]);
    while ~isDone && round < max_rounds
        round = round +1;
        order = zeros(nBlocks ,nBlocks );
        blocks = 1:nBlocks;
        for s = 1:nBlocks
            blocks_s = blocks;
            if s > 1
                blocks_s = blocks_s(~(ismember(blocks_s, order(1:(s-1),1))));
            end
            order(s,1)= blocks_s(randi(numel(blocks_s)));
            success = false; retryCount = 0;
            while ~ success && retryCount < max_retries
                try
                    for i = 2:numel(blocks)
                        blocks_i = blocks;
                        blocks_i = blocks_i(~(ismember(blocks_i,order(s,1:(i-1)))));
                        if s > 1
                            blocks_i = blocks_i(~(ismember(blocks_i,order(1:(s-1),i))));
                        end
                        order(s,i)= blocks_i(randi(numel(blocks_i)));
                    end
                    success = true;
                catch ME
                    retryCount = retryCount +1;
                    %                 disp(['Error occurred: ', ME.message]);
                    %                 disp(['Retrying... Attempt ', num2str(retryCount)]);
                end
            end
            if ~success
                disp('Max retries reached. Exiting...');
                %         else
                %             disp('Operation completed successfully.');
            end
        end
        
        if success
            if round == 1
                all_orders = order;
            else
                all_orders = [all_orders;order];
            end
        end
        % Check Transition Balancing
        transitionMatrix = zeros(nBlocks, nBlocks); % To count how often each transition is used
        for p = 1:size(all_orders, 1)
            trans_Order = all_orders(p, :);
            transitions = [trans_Order(1:end-1)', trans_Order(2:end)']; % Get transitions for the participant
            
            % Update transition frequencies
            for t = 1:size(transitions, 1)
                transitionMatrix(transitions(t, 1), transitions(t, 2)) = ...
                    transitionMatrix(transitions(t, 1), transitions(t, 2)) + 1;
            end
        end
        
        isDone = isequal(transitionMatrix ./ round, perfect_matrix);
    end
    if ~isDone
        disp('Max retries reached. Exiting...');
    else
        disp('Operation completed successfully.');
        save(['all_orders_', sprintf('%02d',file),'.mat'], 'all_orders');
        break
    end
end
end

clear all; 
cd('...\Movie\Block_orders'); 

all01 = load("all_orders_01.mat");
all02 = load("all_orders_02.mat"); 
all03 = load("all_orders_03.mat"); 
all04 = load("all_orders_04.mat"); 
all05 = load("all_orders_05.mat");
all06 = load("all_orders_06.mat"); 
all07 = load("all_orders_07.mat");
all08 = load("all_orders_08.mat");
all09 = load("all_orders_09.mat"); 
all10 = load("all_orders_10.mat"); 
all11 = load("all_orders_11.mat"); 
all12 = load("all_orders_12.mat");
all13 = load("all_orders_13.mat"); 
all14 = load("all_orders_14.mat");
all15 = load("all_orders_15.mat");
all16 = load("all_orders_16.mat"); 
all17 = load("all_orders_17.mat"); 
all18 = load("all_orders_18.mat"); 
all19 = load("all_orders_19.mat");
all20 = load("all_orders_20.mat"); 


all_orders = [all01.all_orders;  all02.all_orders;  all03.all_orders;...
    all04.all_orders;  all05.all_orders;  all06.all_orders;...
    all07.all_orders;  all08.all_orders;  all09.all_orders;...
    all10.all_orders;  all11.all_orders;  all12.all_orders;...
    all13.all_orders;  all14.all_orders;  all15.all_orders;...
    all16.all_orders;  all17.all_orders;  all18.all_orders;...
    all19.all_orders;  all20.all_orders;];
          
all_orders = set_orders;

 transitionMatrix = zeros(6, 6); 
for p = 1:size(all_orders, 1)
    trans_Order = all_orders(p, :);
    transitions = [trans_Order(1:end-1)', trans_Order(2:end)']; % Get transitions for the participant
    
    % Update transition frequencies
    for t = 1:size(transitions, 1)
        transitionMatrix(transitions(t, 1), transitions(t, 2)) = ...
            transitionMatrix(transitions(t, 1), transitions(t, 2)) + 1;
    end
end

disp(transitionMatrix); 

for i = 1:6
    tabulate(all_orders(:,i))
end

    test_equal = zeros(size(all_orders,1),1);
participants = 1:size(all_orders,1);
for i = 1:size(all_orders,1)
    participants_i = participants(~(participants == i));
    for j = participants_i
        test_equal(i) = test_equal(i) + isequal(all_orders(i,:), all_orders(j,:));
    end
end
disp(test_equal'); 

for i = 1:6
    tabulate(all_orders(all_orders(:,1) == 6, i))
end

%% select a random set of 8 options, that are all different and more
% balances in sence of where which block appears

found = false;
max_repeats = 1e4;
repeat = 0; 
while ~found && repeat < max_repeats
    repeat = repeat +1; 
    setOf8 = (randi(20, 8,1)-1)*6+1;
    set_orders = zeros(8*6,6); m = 1; 
    for o = 1:6:(8*6)
        set_orders(o:(o+5),:) = all_orders(setOf8(m):(setOf8(m)+5),:); 
        m = m +1;
    end
    
    
    test_equal = zeros(size(set_orders,1),1);
    participants = 1:size(set_orders,1);
    for i = 1:size(set_orders,1)
        participants_i = participants(~(participants == i));
        for j = participants_i
           test_equal(i) = test_equal(i) + isequal(set_orders(i,:), set_orders(j,:));
        end
    end
    
    howMany = zeros(6,6,6);
    for i = 1:6
          subSet_orders = set_orders(set_orders(:,1) == i, :);
          for j = 1:6
            for k = 1:6
                howMany(k,j,i) = sum(subSet_orders(:,j) == k);
            end
          end
    end
    
    found = ~(any(howMany(:,2:6,:) > 3, 'all')|...
              sum(any(howMany(:,2:6,:) > 2), 'all') > 20|...
              any(test_equal > 0));

end
if ~found
    disp('Max retries reached. Exiting...');
    save('set_orders.mat', 'set_orders'); 
else
    disp('Operation completed successfully.');
end