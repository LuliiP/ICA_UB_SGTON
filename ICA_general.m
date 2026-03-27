%% here we'll write the general code for the ICA 
% new branch for LuisaPenso script 

%% Luisa Cleaning.m
%% path directories

% Luisa Valencia paths
%basePath = 'C:\Users\luisa\OneDrive - Universitat de Barcelona (1)';
%addpath('C:\Users\luisa\OneDrive - Universitat de Barcelona (1)\matlab\eeglab2026.0.0') %to eeglab

% Luisa Penso paths
basePath = '/Users/luisapenso/Library/CloudStorage/OneDrive-UniversitatdeBarcelona/Luisa Maria Valencia Torres''s files - SGTon';
addpath('/Users/luisapenso/Applications/eeglab2025.1.0')
eeglab

SubArray = [1 2 4 5 6 7 8 14];

% Analysis_folder = fullfile(basePath,'SGTon','SCAN EEG');

% Luisa Valencia paths
%Analysis_folder = fullfile(basePath,'SGTon','SCAN EEG');
%ElectrodesTab = fullfile(basePath, 'SGTon', 'SCAN EEG', 'PREP_Reports','Electr.mat');

% Luisa Penso paths
Analysis_folder = fullfile(basePath,'SCAN EEG');
ElectrodesTab = fullfile(basePath, 'SCAN EEG', 'PREP_Reports', 'Electr.mat');

load(ElectrodesTab)

%%
for iSub = 1:length(SubArray)

    subj = sprintf('%02d', SubArray(iSub));

    % Load cleaned continuous rawName = [ subj 'raw68.set'];

    EEG = pop_loadset([subj 'raw68.set'], Analysis_folder);
    orig_chanlocs = EEG.chanlocs;

    % Identify subject row
    subjNum = str2double(subj);
    row = Clean.Subject == subjNum;

    if sum(row)==0
        error('Subject %s not found in Clean table',subj);
    end

    brokenCh = Clean.Broken{row};
    prepLabels = Clean.PREP{row}(:,2);

    % Visual inspection
    EEG_vis = pop_eegfiltnew(EEG,1,45);
    pop_eegplot(EEG_vis,1,1,1);
    fprintf('\nRecommended attention channels:\n')

    if ~isempty(brokenCh)
        fprintf('Broken (hardware): %s\n', strjoin(brokenCh', ', '));
    end

    if ~isempty(prepLabels)
        fprintf('PREP detected: %s\n', strjoin(prepLabels', ', '));
    end

    fprintf('\nThese channels deserve special inspection.\n')
    fprintf('\nEnter channels to remove\n')
    removeCh = input('Channels: ','s'); % write it like this--> VEOGL O2

    % Convert string → cell array of labels
    removeCh = strtrim(split(removeCh, {' ',','}));
    removeCh = removeCh(~cellfun('isempty', removeCh));

    % Overwrite VisualInspection (1-column cell array)
    Clean.VisualInspection{row} = removeCh(:);

    % Remove channels from EEG
    labels = {EEG.chanlocs.labels};
    removeIdx = find(ismember(labels, removeCh));
    EEG = pop_select(EEG,'nochannel',removeIdx);
    EEG = eeg_checkset(EEG);
    
    fprintf('\nBroken (hardware): %s\n', strjoin(brokenCh', ', '));
    %fprintf(fid, '\nBroken (hardware): %s\n', strjoin(brokenCh', ', '));
    %fprintf(fid, '\nPREP detected: %s\n', strjoin(prepLabels', ', '));
    %fprintf(fid, '\nChannels removed (manual): %s\n', strjoin(removeCh', ', '));



    fprintf('Channels removed and table updated\n')


    % save copy file without channels

    cleanName = [ subj 'raw_cleaned.set'];

    pop_saveset(EEG, ...
        'filename',cleanName, ...
        'filepath',Analysis_folder);

    fprintf('Saved clean continuous dataset\n');

    % Prepare ICA copy filtered and rereferenced to average
    % List of possible EOG channel names
    eogCh = {'HEOGL','HEOGR','VEOGU','VEOGL'};

    % Current channel labels
    labels = {EEG.chanlocs.labels};

    % Find which EOG channels are still present
    eogIdx = find(ismember(labels, eogCh));

    % Filter for ICA
    %EEG = pop_firws(EEG, 'fcutoff', 1, 'ftype', 'highpass', 'wtype', 'kaiser', 'warg', 5);
    EEG = pop_firws(EEG, ...
        'fcutoff',1, ...
        'ftype','highpass', ...
        'wtype','kaiser', ...
        'warg',5, ...
        'forder',1650);
    EEG = eeg_checkset(EEG);

    fprintf('Filtering complete\n');

    %EEG = pop_reref(EEG, refChans); %Average reference before ICA
    EEG = pop_reref(EEG, [], 'exclude', eogIdx);

    EEG = eeg_checkset(EEG);

    % COMPUTE RANK


    rankICA = rank(double(eeg_getdatact(EEG)));

    fprintf('Channels: %d | Rank: %d\n',EEG.nbchan,rankICA);
    fprintf(fid, '\nFiltering complete');
    fprintf(fid, '\nRe-referencing data');
    fprintf(fid, '\nChannels: %d | Rank: %d', EEG.nbchan, rankICA);
    fprintf(fid, '\nRunning SOBI ICA...\n');
    % RUN SOBI ICA


    fprintf('Running SOBI ICA...\n');

    EEG = pop_runica(EEG, ...
        'icatype','sobi', ...
        'pca',rankICA);

    EEG = eeg_checkset(EEG);


    % RUN ICLABEL


    EEG = iclabel(EEG);

    EEG = eeg_checkset(EEG);


    nIC = size(EEG.icaweights,1);              % number of ICs actually available
    fprintf('Number of ICs available: %d\n', nIC);

    pop_viewprops(EEG, 0, 1:nIC);              % show only valid ICs
    % MANUAL IC REJECTION

    nIC = size(EEG.icaweights,1);
    fprintf('Number of ICs available: %d\n', nIC);

    badICs_str = input('Enter ICs to reject ([]=none): ','s'); % like this [1 4 5]

    if ~isempty(badICs_str)
        badICs = str2num(badICs_str); %#ok<ST2NM>
        badICs = badICs(badICs >= 1 & badICs <= nIC);   % keep only valid ICs

        EEG = pop_subcomp(EEG, badICs, 0);
        fprintf('Rejected ICs: %s\n', mat2str(badICs));
    else
        badICs = [];
        fprintf('No ICs rejected\n');
    end

    %  SAVE ICA INFO

    ICAinfo.icaweights = EEG.icaweights;
    ICAinfo.icasphere  = EEG.icasphere;
    ICAinfo.icawinv    = EEG.icawinv;
    ICAinfo.badICs     = badICs;

    % Save labels used for ICA
    ICAinfo.chanlabels = {EEG.chanlocs.labels};
    Icaweights = ['subj' subj '_ICAinfo.mat'];

    save(fullfile(ica_folder, ['subj' subj '_ICAinfo.mat']), 'ICAinfo');

    fprintf('ICA weights saved\n');



    % APPLY ICA WEIGHTS TO CLEAN_RAW DATASET

    fprintf('\nApplying ICA weights to full clean_raw dataset...\n');

    cleanName = [ subj 'raw_cleaned.set'];
    EEG_full = pop_loadset('filename', cleanName, 'filepath', Analysis_folder);
    EEG_full = eeg_checkset(EEG_full);

    % Load ICA info
    load(fullfile(ica_folder, ['subj' subj '_ICAinfo.mat']), 'ICAinfo');

    % Select the same channels used during ICA
    EEG_full = pop_select(EEG_full, 'channel', ICAinfo.chanlabels);
    EEG_full = eeg_checkset(EEG_full);

    % Transfer ICA weights
    EEG_full.icaweights = ICAinfo.icaweights;
    EEG_full.icasphere  = ICAinfo.icasphere;
    EEG_full.icawinv    = ICAinfo.icawinv;

    % Reset ICA channel indices
    EEG_full.icachansind = 1:length(ICAinfo.chanlabels);

    fprintf('ICA weights applied successfully.\n');


    % REMOVE BAD ICS FROM FULL DATA

    if ~isempty(ICAinfo.badICs)

        fprintf('Removing ICs from continuous dataset...\n');

        EEG_full = pop_subcomp(EEG_full, ICAinfo.badICs, 0);

    end


    % INTERPOLATE REMOVED CHANNELS

    fprintf('Interpolating removed channels...\n');

    % Load original montage
    %orig = pop_loadset('filename', fileName, 'filepath', Raquel_folder);

    % IMPORTANT: clear ICA fields before interpolation
    EEG_full.icaweights = [];
    EEG_full.icasphere  = [];
    EEG_full.icawinv    = [];
    EEG_full.icachansind = [];
    EEG_full.icaact     = [];

    scalp_idx = find(strcmp({orig_chanlocs.type}, 'EEG'));
    EEG_full = pop_interp(EEG_full, orig_chanlocs(scalp_idx), 'spherical');%EEG_full = pop_interp(EEG_full, orig_chanlocs, 'spherical');

    EEG_full = eeg_checkset(EEG_full);

    fprintf('Missing channels interpolated.\n');


    % SAVE FINAL file
    fileName = [subj '_ALLCOND_trimmed.set'];
    rawName = [ subj 'raw68.set'];
    finalName = [subj '_ICAclean_full.set'];



    % Build the comment block
    myCommentBlock = [ ...
        '1. All original CDT files for this subject were merged into a single continuous dataset.' 10 ...
        '   Source files located in:' 10 ...
        % '      ' fullfile('C:\Users\luisa\OneDrive - Universitat de Barcelona (1)\SGTon\Raw_Data\', ['S' subj '_C01']) 10 ...
        % '      ...' 10 ...
        % '      ' fullfile('C:\Users\luisa\OneDrive - Universitat de Barcelona (1)\SGTon\Raw_Data\', ['S' subj '_C02']) 10 ...

        % NEW
        '      ' fullfile(basePath, 'Raw_Data', ['S' subj '_C01']) 10 ...
        '      ' fullfile(basePath, 'Raw_Data', ['S' subj '_C02']) 10 ...

        '   Output: ' fileName 10 10 ...
        '2. Manual artifact segments were removed based on predefined criteria.' 10 10 ...
        '3. PREP pipeline was executed in detection-only mode (no interpolation, no RANSAC).' 10 ...
        '   Output: ' rawName 10 10 ...
        '4. Channels flagged by hardware notes, PREP, and visual inspection were removed.' 10 ...
        '   Output: ' cleanName 10 10 ...
        '5. A filtered (1 Hz) and average-referenced copy was created for ICA computation.' 10 ...
        '   SOBI ICA was computed on this filtered copy.' 10 ...
        '   ICA weights, sphere, and inverse weights were saved to: ' Icaweights 10 ...
        '   These ICA weights were then applied back to the unfiltered clean dataset.' 10 10 ...
        '6. After ICA component rejection, removed channels were interpolated using the original montage.' 10 ...
        '   Final output: ' finalName 10 ...
        ];

    % Append comments correctly
    EEG_full.comments = pop_comments(EEG_full.comments, '', myCommentBlock, 1);

    % IMPORTANT TO RUN THESE LINES (if error in myCommentBlock, comment it out with %)
    pop_saveset(EEG_full, ...
        'filename', finalName, ...
        'filepath', Analysis_folder);

    fprintf('Saved final ICA-cleaned dataset: %s\n', finalName);

end

% check the file !!



%%
basePath = 'C:\Users\luisa\OneDrive - Universitat de Barcelona (1)';
addpath('C:\Users\luisa\OneDrive - Universitat de Barcelona (1)\matlab\eeglab2026.0.0') %to eeglab
eeglab
close all
Analysis_folder = fullfile(basePath,'SGTon','SCAN EEG');
ica_folder = fullfile(basePath, 'SGTon', 'SCAN EEG', 'ICA_w');
if ~exist(ica_folder, 'dir'), mkdir(ica_folder); end
ElectrodesTab = fullfile(basePath, 'SGTon', 'SCAN EEG', 'PREP_Reports','Electr.mat');

logFile = fullfile(ica_folder,'ICA_log.txt');
fid = fopen(logFile,'a');

fprintf(fid,'\n===================================');
fprintf(fid,'\n%s',datestr(now));
fprintf(fid,'\nPREPROCESSING PIPELINE\n');

%load(ElectrodesTab)
SubArray = [1 2 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 21 22 24 25 27 28 30 32 33 34 35] %3 missing
%%P5
for iSub = 1:length(SubArray)
  load(ElectrodesTab)

    subj = sprintf('%02d', SubArray(iSub));

    % Load cleaned continuous 
    rawName = [ subj 'raw68.set'];

    EEG = pop_loadset('filename', rawName, 'filepath', Analysis_folder); %EEG = pop_loadset([rawName, Analysis_folder);
    orig_chanlocs = EEG.chanlocs;

    % Identify subject row
    subjNum = str2double(subj);
    row = Clean.Subject == subjNum;

    if sum(row)==0
        error('Subject %s not found in Clean table',subj);
    end

    brokenCh = Clean.Broken{row};
    prepLabels = Clean.PREP{row}(:,2);

    % Visual inspection
    EEG_vis = pop_eegfiltnew(EEG,1,45);
    pop_eegplot(EEG_vis,1,1,1);
    fprintf('\nRecommended attention channels:\n')

    if ~isempty(brokenCh)
        fprintf('Broken (hardware): %s\n', strjoin(brokenCh', ', '));
    end

    if ~isempty(prepLabels)
        fprintf('PREP detected: %s\n', strjoin(prepLabels', ', '));
    end

    fprintf('\nThese channels deserve special inspection.\n')
    fprintf('\nEnter channels to remove\n')
    removeCh = input('Channels: ','s'); % write it like this--> VEOGL O2

    % Convert string → cell array of labels
    removeCh = strtrim(split(removeCh, {' ',','}));
    removeCh = removeCh(~cellfun('isempty', removeCh));

    % Overwrite VisualInspection (1-column cell array)
    Clean.VisualInspection{row} = removeCh(:);

    % Remove channels from EEG
    labels = {EEG.chanlocs.labels};
    removeIdx = find(ismember(labels, removeCh));
    EEG = pop_select(EEG,'nochannel',removeIdx);
    EEG = eeg_checkset(EEG);
    
    fprintf('\nBroken (hardware): %s\n', strjoin(brokenCh', ', '));
    %fprintf(fid, '\nBroken (hardware): %s\n', strjoin(brokenCh', ', '));
    %fprintf(fid, '\nPREP detected: %s\n', strjoin(prepLabels', ', '));
    %fprintf(fid, '\nChannels removed (manual): %s\n', strjoin(removeCh', ', '));



    fprintf('Channels removed and table updated\n')


    % save copy file without channels

    cleanName = [ subj 'raw_cleaned.set'];

    pop_saveset(EEG, ...
        'filename',cleanName, ...
        'filepath',Analysis_folder);

    fprintf('Saved clean continuous dataset\n');

    % Prepare ICA copy filtered and rereferenced to average
    % List of possible EOG channel names
    eogCh = {'HEOGL','HEOGR','VEOGU','VEOGL'};

    % Current channel labels
    labels = {EEG.chanlocs.labels};

    % Find which EOG channels are still present
    eogIdx = find(ismember(labels, eogCh));

    % Filter for ICA
    %EEG = pop_firws(EEG, 'fcutoff', 1, 'ftype', 'highpass', 'wtype', 'kaiser', 'warg', 5);
    EEG = pop_firws(EEG, ...
        'fcutoff',1, ...
        'ftype','highpass', ...
        'wtype','kaiser', ...
        'warg',5, ...
        'forder',1650);
    EEG = eeg_checkset(EEG);

    fprintf('Filtering complete\n');

    %EEG = pop_reref(EEG, refChans); %Average reference before ICA
    EEG = pop_reref(EEG, [], 'exclude', eogIdx);

    EEG = eeg_checkset(EEG);

    % COMPUTE RANK


    rankICA = rank(double(eeg_getdatact(EEG)));

    fprintf('Channels: %d | Rank: %d\n',EEG.nbchan,rankICA);
    fprintf(fid, '\nFiltering complete');
    fprintf(fid, '\nRe-referencing data');
    fprintf(fid, '\nChannels: %d | Rank: %d', EEG.nbchan, rankICA);
    fprintf(fid, '\nRunning SOBI ICA...\n');
    % RUN SOBI ICA


    fprintf('Running SOBI ICA...\n');

    EEG = pop_runica(EEG, ...
        'icatype','sobi', ...
        'pca',rankICA);

    EEG = eeg_checkset(EEG);


    % RUN ICLABEL

EEG = pop_iclabel(EEG, 'default');
    ICclass  = EEG.etc.ic_classification.ICLabel.classifications;
    IClabels = EEG.etc.ic_classification.ICLabel.classes;
    
    idxEye = find(strcmpi(IClabels, 'Eye'));
    idxMus = find(strcmpi(IClabels, 'Muscle'));
    idxBrn = find(strcmpi(IClabels, 'Brain'));

    % Umbrales para sugerencias
    thr_eye    = 0.60;
    thr_muscle = 0.80;
    
    candidateICs = find( ...
        (ICclass(:,idxEye) >= thr_eye | ...
         ICclass(:,idxMus) >= thr_muscle) & ...
         ICclass(:,idxBrn) < 0.50 );

    fprintf('\n--- Automatic Suggestions ---\n');
    fprintf('Eye >= %.2f | Muscle >= %.2f | Brain < 0.50\n', thr_eye, thr_muscle);
    fprintf('Suggested ICs to reject: %s\n', mat2str(candidateICs'));

    % Visualization 
    pop_eegplot(EEG, 0, 1, 1); 
    fprintf('Look the different componets overtime...\n');

    % Manual input
    nIC = size(EEG.icaweights, 1);
    badICs_str = input('Enter ICs to reject (like  [1 4 5] o []): ', 's');
    
    if ~isempty(badICs_str)
        badICs = str2num(badICs_str); %#ok<ST2NM>
        badICs = badICs(badICs >= 1 & badICs <= nIC); % Validación de rango
        
        % Aplicar rechazo
        EEG = pop_subcomp(EEG, badICs, 0);
        fprintf('Rejected ICs: %s\n', mat2str(badICs));
    else
        badICs = [];
        fprintf('No ICs rejected\n');
    end


    %  SAVE ICA INFO

    ICAinfo.icaweights = EEG.icaweights;
    ICAinfo.icasphere  = EEG.icasphere;
    ICAinfo.icawinv    = EEG.icawinv;
    ICAinfo.badICs     = badICs;

    % Save labels used for ICA
    ICAinfo.chanlabels = {EEG.chanlocs.labels};
    Icaweights = ['subj' subj '_ICAinfo.mat'];

    save(fullfile(ica_folder, ['subj' subj '_ICAinfo.mat']), 'ICAinfo');

    fprintf('ICA weights saved\n');



    % APPLY ICA WEIGHTS TO CLEAN_RAW DATASET

    fprintf('\nApplying ICA weights to full clean_raw dataset...\n');

    cleanName = [ subj 'raw_cleaned.set'];
    EEG_full = pop_loadset('filename', cleanName, 'filepath', Analysis_folder);
    EEG_full = eeg_checkset(EEG_full);

    % Load ICA info
    load(fullfile(ica_folder, ['subj' subj '_ICAinfo.mat']), 'ICAinfo');

    % Select the same channels used during ICA
    EEG_full = pop_select(EEG_full, 'channel', ICAinfo.chanlabels);
    EEG_full = eeg_checkset(EEG_full);

    % Transfer ICA weights
    EEG_full.icaweights = ICAinfo.icaweights;
    EEG_full.icasphere  = ICAinfo.icasphere;
    EEG_full.icawinv    = ICAinfo.icawinv;

    % Reset ICA channel indices
    EEG_full.icachansind = 1:length(ICAinfo.chanlabels);

    fprintf('ICA weights applied successfully.\n');


    % REMOVE BAD ICS FROM FULL DATA

    if ~isempty(ICAinfo.badICs)

        fprintf('Removing ICs from continuous dataset...\n');

        EEG_full = pop_subcomp(EEG_full, ICAinfo.badICs, 0);

    end


    % INTERPOLATE REMOVED CHANNELS

    fprintf('Interpolating removed channels...\n');

    % Load original montage
    %orig = pop_loadset('filename', fileName, 'filepath', Raquel_folder);

    % IMPORTANT: clear ICA fields before interpolation
    EEG_full.icaweights = [];
    EEG_full.icasphere  = [];
    EEG_full.icawinv    = [];
    EEG_full.icachansind = [];
    EEG_full.icaact     = [];

    scalp_idx = find(strcmp({orig_chanlocs.type}, 'EEG'));
    EEG_full = pop_interp(EEG_full, orig_chanlocs(scalp_idx), 'spherical');%EEG_full = pop_interp(EEG_full, orig_chanlocs, 'spherical');

    EEG_full = eeg_checkset(EEG_full);

    fprintf('Missing channels interpolated.\n');


    % SAVE FINAL file
    fileName = [subj '_ALLCOND_trimmed.set'];
    finalName = [subj '_ICAclean_full.set'];
    save(ElectrodesTab, 'Clean');
    fprintf('Table "Clean" updated and saved for subject %s\n', subj);


    % Build the comment block
    myCommentBlock = [ ...
        '1. All original CDT files for this subject were merged into a single continuous dataset.' 10 ...
        '   Source files located in:' 10 ...
        % '      ' fullfile('C:\Users\luisa\OneDrive - Universitat de Barcelona (1)\SGTon\Raw_Data\', ['S' subj '_C01']) 10 ...
        % '      ...' 10 ...
        % '      ' fullfile('C:\Users\luisa\OneDrive - Universitat de Barcelona (1)\SGTon\Raw_Data\', ['S' subj '_C02']) 10 ...

        % NEW
        '      ' fullfile(basePath, 'Raw_Data', ['S' subj '_C01']) 10 ...
        '      ' fullfile(basePath, 'Raw_Data', ['S' subj '_C02']) 10 ...

        '   Output: ' fileName 10 10 ...
        '2. Manual artifact segments were removed based on predefined criteria.' 10 10 ...
        '3. PREP pipeline was executed in detection-only mode (no interpolation, no RANSAC).' 10 ...
        '   Output: ' rawName 10 10 ...
        '4. Channels flagged by hardware notes, PREP, and visual inspection were removed.' 10 ...
        '   Output: ' cleanName 10 10 ...
        '5. A filtered (1 Hz) and average-referenced copy was created for ICA computation.' 10 ...
        '   SOBI ICA was computed on this filtered copy.' 10 ...
        '   ICA weights, sphere, and inverse weights were saved to: ' Icaweights 10 ...
        '   These ICA weights were then applied back to the unfiltered clean dataset.' 10 10 ...
        '6. After ICA component rejection, removed channels were interpolated using the original montage.' 10 ...
        '   Final output: ' finalName 10 ...
        ];

    % Append comments correctly
    EEG_full.comments = pop_comments(EEG_full.comments, '', myCommentBlock, 1);

    % IMPORTANT TO RUN THESE LINES (if error in myCommentBlock, comment it out with %)
    pop_saveset(EEG_full, ...
        'filename', finalName, ...
        'filepath', Analysis_folder);

    fprintf('Saved final ICA-cleaned dataset: %s\n', finalName);

end
fclose(fid)

% check the file !!

