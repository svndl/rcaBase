function rca_struct = rcaSweep(path_names, rca_bins, rca_freqs, rca_conds, rca_trials, rca_subs, nReg, nComp, dataType, chanToCompare, forceReload)
% perform RCA on sweep SSVEP data exported to RLS or DFT format
%
% [rcaData,W,A,noiseData,compareData,compareNoiseData,freqIndices,binIndices]=RCASWEEP(path_names,[rca_bins],[rca_freqs],[rca_conds],[rca_trials],[NREG],[NCOMP],[DATATYPE],[COMPARECHAN],[SHOW],[RCPLOTSTYLE])
%
% INPUTS:
% path_names:   cell vector of string directory names housed in DFT_c00x.txt or RLS_c00x.txt exports
%                   for example, path_names={'/Volumes/Denali_4D2/rca/s001/','/Volumes/Denali_4D2/rca/s002/','/Volumes/Denali_4D2/rca/s003/'}
% rca_bins:    vector of bin indices indices to to use for rca training
%                   [ defaults to bin 0 or average across bins ] 
% rca_freqs:   vector of frequency indices to to use for rca training
%                   [ defaults to all available harmonics ] 
% rca_conds:   vector of conditions to use for rca training
% rca_trials:  vector of indices of trials to use for rca training
% rca_subs:    vector of subjects to use for rca training
% nReg:        RCA regularization parameter (defaults to 9)
% nComp:       number of RCs to retain (defaults to 3)
% dataType:    can be 'DFT' or 'RLS'
% compareChan: comparison channel index between 1 and the total number
%              of channels in the specified dataset [75]
% forceReload: true/[false], if true, reload the data text files and generate
%              new .mat files for quick loading of the data
%              set to true if you have re-exported the data
%
% OUTPUTS:
% rcaData: returned by rcaRun
% W: linear transformation matrix to go from sensor space to RC-space
% A: linear transformation matrix to go from RC-space space to sensor space
% covData: the covariance data used to learn the RCs
% noiseData: struct containing fields for the rca space representation of frequency side bands
% compareData: rcaData for the comparison channel only (if supplied)
% compareNoiseData: noiseData for the comparison channel only (if supplied)
% rcaSettings: struct containing settings used to run RCA & indices
%              necessary for unpacking rcaData later on
%
% rcaData is a (NumberConditions x NumberSubjects) cell array. Each cell
% contains a (length(rca_freqs)*length(rca_bins)*2 x nComp x NumberTrialsPerSubject)
% matrix containing the real and imaginary components of each RC.
%
% Jacek P. Dmochowski, 2015, report bugs to dmochowski@gmail.com
% Edited by HEG 07/2015

if nargin<11 || isempty(forceReload); forceReload = false; end
if nargin<10 || isempty(chanToCompare); chanToCompare = 75; end
if nargin<9 || isempty(dataType), dataType = 'RLS'; end
if nargin<8 || isempty(nComp), nComp=3; end
if nargin<7 || isempty(nReg), nReg=9; end
if nargin<6 || isempty(rca_subs), rca_subs = []; end
if nargin<5 || isempty(rca_trials), rca_trials = []; end
if nargin<4 || isempty(rca_conds), rca_conds=[]; end
if nargin<3 || isempty(rca_freqs), rca_freqs=[]; end
if nargin<2 || isempty(rca_bins), rca_bins=[]; end
if nargin<1, error('Must specify at least one input argument'); end

%% if pathanmes is a string, convert to cell
if ~iscell(path_names)
    try
        path_names=mat2cell(path_names);
    catch
        error('unable to parse path_names: check that it is a cell array of strings');
    end
end
if isempty(path_names)
    error('Variable path_names is empty.');
end

%% read in signal and noise data
nSubjects=numel(path_names);
all_data={};
cellNoiseData1={};
cellNoiseData2={};
subFreqIdx=cell(nSubjects,1);
subBinIdx=cell(nSubjects,1);
subFreqLabels=cell(nSubjects,1);
subBinLabels=cell(nSubjects,1);
fprintf('Reading in sensor data from provided path names...\n')
s = 1;
assigned = false(1,length(rca_conds));
while (s <= nSubjects)
    sourceDataFileName = sprintf('%s/sourceData_%s.mat',path_names{s},dataType);
    if isempty(dir(sourceDataFileName)) || forceReload
        createSourceDataMat(path_names{s});
    end
    [signalData,noise1,noise2,subFreqIdx{s},subBinIdx{s},subFreqLabels{s},subBinLabels{s}] = selectDataForTraining(sourceDataFileName);
    all_data(:,s)=signalData;
    cellNoiseData1(:,s)=noise1;
    cellNoiseData2(:,s)=noise2;
    train_data(:,s) = selectDataForTraining(sourceDataFileName,rca_bins,rca_freqs,rca_conds,rca_trials);
    
    for c = 1:length(rca_conds)
        if ~any(isnan(subFreqIdx{s}{c}))
            if assigned(c)
                % test against master list
                isMatched  = [isequal(freqIndices{c},subFreqIdx{s}{c}), ...
                              isequal(binIndices{c},subBinIdx{s}{c}), ...
                              isequal(freqLabels{c},subFreqLabels{s}{c}), ...
                              isequal(binLabels{c},subBinLabels{s}{c})];
                if sum( ~isMatched ) > 0
                    error('Frequency and bin indices vary across subjects: check consistency of DFT/RLS exports\n.');
                else
                    clear isMatched;
                end
            else
                assigned(c) = true;
                % assign to master list
                freqIndices{c} = subFreqIdx{s}{c};
                binIndices{c} = subBinIdx{s}{c};
                freqLabels{c} = subFreqLabels{s}{c};
                binLabels{c} = subBinLabels{s}{c};
            end
        else
        end
    end
    fprintf('Done selecting data for subject %d/%d: %s\n',s,nSubjects,path_names{s});
    s = s + 1;
end
% now check that freqs and bins are consistent across conditions
isConsistent(1) = all(cellfun(@(x) isequal(freqIndices{1},x), freqIndices));
isConsistent(2) = all(cellfun(@(x) isequal(binIndices{1},x), binIndices));
isConsistent(3) = all(cellfun(@(x) isequal(freqLabels{1},x), freqLabels));
isConsistent(4) = all(cellfun(@(x) isequal(binLabels{1},x), binLabels));

if any(~isConsistent)
    msg = 'Frequency and bin indices and/or labels vary across conditions \n';
    error(msg);
else
    freqIndices = freqIndices{1};
    binIndices = binIndices{1};
    freqLabels = freqLabels{1};
    binLabels = binLabels{1};
end
nChannels = size(all_data{1,1},2);

if chanToCompare > nChannels
    error_msg = sprintf('\n The requested comparision channel (%d) exceeds the maximum number of channels in the dataset (%d), ignoring request...\n', chanToCompare, nChannels);
    error(error_msg);
else
end

%% DEAL WITH BINS
if any(contains(binLabels, 'NaN'))
    nan_idx = contains(binLabels, 'NaN') == 1;
    binLabels(nan_idx) = arrayfun(@(x) sprintf('bin%02d', (x)), 1:sum(nan_idx), 'uni' , false);
else
end

if any(binIndices == 0)
    binIndices = binIndices + 1; % turn into actual indices
    % move averages to the end of the stack
    for b = 1:length(binIndices)
        if binIndices(b) == 11
            new_idx(b) = b - 10;
        else
            new_idx(b) = b + 1;
        end
    end
    new_idx = [new_idx, new_idx+length(new_idx)];
    all_data = cellfun(@(x) x(new_idx, :, :), all_data, 'uni', false);
    % note: not necessary to reorder train_data, as this 
    % does not influence the output
    binLabels = binLabels(new_idx(1:length(binLabels)));
    % make sure rca bins = 0 refer to the now 11th index ("ave")
    rca_bins(rca_bins == 0) = 11;
    clear new_idx;
else
end

%% DEAL WITH EMPTY RCA INPUTS
if isempty(rca_conds)
    rca_conds = 'all';
else
end
if isempty(rca_trials)
    rca_trials = 'all';
else
end 
if isempty(rca_bins)
    rca_bins = 'all';
else
end
if isempty(rca_freqs)
    rca_freqs = 'all';
else
end

% make new array that only has training subjects

if ~isempty(rca_subs)
    if any(~ismember(rca_subs, 1:size(train_data, 2)))
        msg = '\n subject index must not be larger than number of subjects \n';
        error(msg)
    else
    end
    train_data = train_data(:,rca_subs);
else
    rca_subs = 'all';
end

%% run RCA
fprintf('Running RCA...\n');
warning('off','all')
[~, W, A, Rxx, Ryy, Rxy, dGen] = rcaRun(train_data, nReg, nComp); 
covData.Rxx = Rxx; 
covData.Ryy = Ryy;
covData.Rxy = Rxy;
covData.sortedGeneralizedEigenValues = dGen;
% rca was trained on 'train_data', 
% but now project all the data through those components
rcaData = rcaProject(all_data, W);
noiseData.lowerSideBand=rcaProject(cellNoiseData1, W); 
noiseData.higherSideBand=rcaProject(cellNoiseData2, W);

%% CREATE A COMPONENT OF JUST ONE CHANNEL, FOR PERFORMANCE EVALUATION
wComparison=zeros(nChannels,1); wComparison(chanToCompare)=1; 
comparisonData=rcaProject(all_data, wComparison); 
comparisonNoiseData.lowerSideBand =rcaProject(cellNoiseData1,wComparison); 
comparisonNoiseData.higherSideBand =rcaProject(cellNoiseData2,wComparison); 
warning('on','all')

%% package some of the necessary variables for plotting/grouping data later on
rcaSettings.rcaBins = rca_bins;
rcaSettings.rcaFreqs = rca_freqs;
rcaSettings.rcaConds = rca_conds;
rcaSettings.rcaTrials = rca_trials;
rcaSettings.rcaSubs = rca_subs;
rcaSettings.nReg = nReg;
rcaSettings.nComp = nComp;
rcaSettings.chanToCompare = chanToCompare;
rcaSettings.dataFolders = path_names;
rcaSettings.freqLabels = freqLabels; 
rcaSettings.binLabels = binLabels;
rcaSettings.dataType = dataType;
runDate = datestr(clock,26);
runDate(strfind(runDate,'/')) ='';
rcaSettings.runDate = runDate; % on what date was this RCA run?

%% generate final output struct array, organized by frequency in input data

%% THIS SECTION NEEDS TO BE MODERATED
all_freqs = unique(freqIndices);

for f = 1:length(all_freqs)
    rca_struct(f).W = W;
    rca_struct(f).A = A;
    rca_struct(f).covData = covData;
    fIdx = repmat(freqIndices==all_freqs(f),2,1); % repmat because the first half is real, second half is imag with same ordering
    rca_struct(f).data = cellfun(@(x) x(fIdx,:,:),rcaData,'uni',false);
    rca_struct(f).noiseData.lowerSideBand = cellfun(@(x) x(fIdx,:,:),noiseData.lowerSideBand,'uni',false);
    rca_struct(f).noiseData.higherSideBand = cellfun(@(x) x(fIdx,:,:),noiseData.higherSideBand,'uni',false);
    rca_struct(f).comparisonData = cellfun(@(x) x(fIdx,:,:),comparisonData,'uni',false);
    rca_struct(f).comparisonNoiseData.lowerSideBand = cellfun(@(x) x(fIdx,:,:),comparisonNoiseData.lowerSideBand,'uni',false);
    rca_struct(f).comparisonNoiseData.higherSideBand = cellfun(@(x) x(fIdx,:,:),comparisonNoiseData.higherSideBand,'uni',false);
    rca_struct(f).inputData = cellfun(@(x) x(fIdx,:,:),all_data,'uni',false);
    % output indices that match the data
    rcaSettings.freqIndices = freqIndices(freqIndices==all_freqs(f));
    rcaSettings.binIndices = binIndices(freqIndices==all_freqs(f));
    % put into output struct
    rca_struct(f).settings = orderfields(rcaSettings);
    % fields for putting aggregate data later (see aggregateData fxn)
    rca_struct(f).mean = struct([]);
    rca_struct(f).subjects = struct([]); 
    rca_struct(f).stats = struct([]); 
end