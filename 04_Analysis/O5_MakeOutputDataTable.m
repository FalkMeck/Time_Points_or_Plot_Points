%%, CREATE RESULTS TABLE Á LA SPM12 FROM AFNI REPORTS
study_dir = '...\NIFTI\_AFNI_Analysis';
cd(study_dir);

contrasts = {'01_Shot_Scene', '02_04s_12s', '03_04s_36s', '04_12s_36s', ...
    'INT01_12sSceneShot_04sSceneShot','INT02_36sSceneShot_04sSceneShot','INT03_36sSceneShot_12sSceneShot'};
for con = 1:numel(contrasts)
    
    %% Read Files
    reportDir = fullfile(study_dir, contrasts{con}, 'ClusterizeOutput_report');
    report = fullfile(reportDir, [contrasts{con},'_report.1D']);
    
    fid = fopen(report, 'r');
    
    headerLines = {};
    dataLines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(strtrim(line), '#')
            headerLines{end+1} = line;  % Collect header
        elseif ~isempty(strtrim(line))
            dataLines{end+1} = line;    % Collect data
        end
    end
    fclose(fid);
    
    %% Ectract Meta data and clean it up a bit
    metadata = struct();
    for i = 1:length(headerLines)
        line = strtrim(headerLines{i});
        % Match lines like: #[ Key = Value ]
        if startsWith(line, '#[')
            tokens = regexp(line, '#\[\s*(.*?)\s*=\s*(.*?)\s*\]', 'tokens');
            if ~isempty(tokens)
                key = matlab.lang.makeValidName(strtrim(tokens{1}{1}));
                value = strtrim(tokens{1}{2});
                metadata.(key) = value;
            end
        end
    end
    
    cleanMeta = struct();
    
    fields = fieldnames(metadata);
    for i = 1:numel(fields)
        key = fields{i};
        val = strtrim(metadata.(key));
        
        % Remove surrounding brackets or quotes
        val = regexprep(val, "^[\[\']+|[\]']+$", "");
        
        % Try to convert to number if it's numeric
        numVal = str2double(val);
        if ~isnan(numVal)
            cleanMeta.(key) = numVal;
        elseif contains(val, ',')  % might be a list
            cleanMeta.(key) = strsplit(val, ',');
        else
            cleanMeta.(key) = val;  % leave as string
        end
    end
    
    save([reportDir,filesep,'MetaData_',contrasts{con},'.mat'],'cleanMeta');
    
    %% Extract header/var names for data table
    headerCols = {};
    for i = 1:length(headerLines)
        line = strtrim(headerLines{i});
        if contains(line, 'Volume') && contains(line, 'CM LR')
            headerCols = strsplit(strtrim(erase(line, '#')));
            break;
        end
    end
    headerNames = {headerCols{1}, [headerCols{2},'_',headerCols{3}],...
        [headerCols{4},'_',headerCols{5}], [headerCols{6},'_',headerCols{7}],...
        headerCols{8:15}, [headerCols{16},'_',headerCols{17}], [headerCols{18},'_',headerCols{19}],...
        [headerCols{20},'_',headerCols{21}], [headerCols{22},'_',headerCols{23}]};
    
    %% make data table inclusing local maxima
    numData = cell2mat(cellfun(@(l) sscanf(l, '%f')', dataLines, 'UniformOutput', false));
    data = reshape(numData, numel(headerNames), numel(dataLines))';
    
    reportCont = dir(reportDir);
    reportNames = {reportCont.name};
    clustPeakNames = reportNames(endsWith(reportNames, 'peaks_coords.txt'));
    
    if numel(clustPeakNames) == numel(dataLines)
        dataPeaks = [];
        for i = 1:numel(clustPeakNames)
            %disp(i);
            dataPeaks = [dataPeaks; data(i,:)];
            peakName = ['Con_', contrasts{con},'_cluster_', num2str(i),'_peaks_coords.txt'];
            filename = fullfile(reportDir, peakName);
            fid = fopen(filename, 'r');
            
            if fid == -1
                error('Could not open file: %s', filename);
            end
            
            rawLines = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            rawLines = rawLines{1};
            
            if isempty(rawLines)
                localPeaks = NaN;  % File is empty
            else
                % Split each line into numbers
                localPeaks = cellfun(@(line) str2double(strsplit(strtrim(line))), rawLines, 'UniformOutput', false);
                localPeaks = vertcat(localPeaks{:});
            end
            
            
            if ~isnan(localPeaks)
                localAdd = NaN(size(localPeaks,1), size(data,2));
                localAdd(:, ismember(headerNames,'Max_Int')) = localPeaks(:,4);
                localAdd(:, ismember(headerNames,'MI_LR')|...
                    ismember(headerNames,'MI_PA')|...
                    ismember(headerNames,'MI_IS')) = localPeaks(:,1:3);
                
                [~, idx] = sort(abs(localAdd(:, ismember(headerNames,'Max_Int'))),'descend');
                localAdd = localAdd(idx, :);
                
                dataPeaks = [dataPeaks; localAdd];
            end
        end
    else
        disp('error msg');
    end
    dataTable = array2table(dataPeaks,'VariableNames', headerNames);
    
    % Make a copy of the table for export
    exportTable = dataTable;
    
    % Convert all numeric columns' NaNs to empty strings
    for i = 1:size(exportTable,2)
        col = exportTable.(i);
        if isnumeric(col)
            % Convert to cell array with empty strings for NaNs
            colCell = num2cell(col);
            colCell(isnan(col)) = {''};
            exportTable.(i) = colCell;  % now it's a cell array
        end
    end
    outfile = fullfile(reportDir, ['ClustTab_',contrasts{con},'.xlsx']);
    % Write to Excel
    writetable(exportTable, outfile);
    
    
end
