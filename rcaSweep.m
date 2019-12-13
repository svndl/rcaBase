function rca_struct = rcaSweep(path_names, rca_bins, rca_harms, rca_conds, rca_trials, rca_subs, n_reg, n_comp, data_type, compare_ch, force_reload, return_all)
% perform RCA on sweep SSVEP data exported to RLS or DFT format
%
% rca_struct = rcaSweep(path_names, [rca_bins, rca_harms, rca_conds, rca_trials, rca_subs, n_reg, n_comp, data_type, compare_ch, force_reload])
%
% INPUTS:
%   path_names: cell vector of string directory names housed in DFT_c00x.txt or RLS_c00x.txt exports
%                   for example, path_names={'/Volumes/Denali_DATA1/rca/s001/','/Volumes/Denali_DATA1/rca/s002/','/Volumes/Denali_DATA1/rca/s003/'}
%   rca_bins:   vector of indices of bins to use for rca training
%                   [ defaults to 0 which is the average across bins ] 
%   rca_harms:  vector of indices of harmonics to use for rca training
%                   [ defaults to all harmonics in the data set ] 
%   rca_conds:  vector of conditions to use for rca training
%                   [ defaults to all conditions in the data set ] 
%   rca_trials: vector of indices of trials to use for rca training
%                   [ defaults to all trials in the data set ]
%   rca_subs:   vector of subjects to use for rca training
%                   [ defaults to all subjects in the data set ]
%   n_reg:      RCA regularization parameter [9]
%
%   n_comp:     number of RCs to retain [3]
%
%   data_type:  can be 'DFT' or ['RLS']
%
%   compare_ch: comparison channel index between 1 and the total number
%                   of channels in the specified dataset [75]
%   force_reload: true/[false], if true, reload the data text files and generate
%                   new .mat files for quick loading of the data
%                   set to true if you have re-exported the data
%
%   return_all: [true]/false, if true return all of the data in the dataset,
%                   if false, return only data used to train rca
% OUTPUTS: 
%   rca_struct: length(rca_harms) x 1 array of structs, one for each harmonic, 
%   with the following subfields:
%   - rca_data: a conditions x subjects cell array. Each cell contains a 
%           (n_bins+1)*2 x n_comp+1 x trials matrix
%           containing the real and imaginary components of each RC.
%           The comparison data is added at the n_comp+1 position in the matrix
%           (note: by default this is all the data in the dataset, not just the data used for RCA)
%
%   - noiseLower and noiseHigher:
%           two fields containing matrices same shape as data, 
%           with the "noise" in the lower and higher frequency side bands 
%           passed through the components and comparison channel
%
%   - input_data: the input data 
%           (note: by default this is all the data in the dataset, not just the data used for RCA)
%
%   - W: linear transformation matrix to go from sensor space to RC-space
%
%   - A: linear transformation matrix to go from RC-space space to sensor space
%
%   - COV: the covariance data used to learn the RCs
%
%   - settings: struct containing settings used to run RCA and
%              indices necessary for unpacking data later on
%
% created by JP Dmochowski and H Gerhard
% updated and maintained by PJ Kohler [pjkohl3r@gmail.com]

if nargin<12 || isempty(return_all); return_all = true; end
if nargin<11 || isempty(force_reload); force_reload = false; end
if nargin<10 || isempty(compare_ch); compare_ch = 75; end
if nargin<9 || isempty(data_type), data_type = 'RLS'; end
if nargin<8 || isempty(n_comp), n_comp=3; end
if nargin<7 || isempty(n_reg), n_reg=9; end
if nargin<6 || isempty(rca_subs), rca_subs = []; end
if nargin<5 || isempty(rca_trials), rca_trials = []; end
if nargin<4 || isempty(rca_conds), rca_conds = []; end
if nargin<3 || isempty(rca_harms), rca_harms = []; end
if nargin<2 || isempty(rca_bins), rca_bins = 0; end
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
    sourceDataFileName = sprintf('%s/sourceData_%s.mat',path_names{s},data_type);
    if isempty(dir(sourceDataFileName)) || force_reload
        createSourceDataMat(path_names{s});
    end
    if return_all
        [all_data(:,s), cellNoiseData1(:,s), cellNoiseData2(:,s), subFreqIdx{s}, subBinIdx{s}, subFreqLabels{s}, subBinLabels{s}] = selectDataForTraining(sourceDataFileName);
        train_data(:,s) = selectDataForTraining(sourceDataFileName, rca_bins,rca_harms,rca_conds,rca_trials);
    else
        [all_data(:,s), cellNoiseData1(:,s), cellNoiseData2(:,s), subFreqIdx{s}, subBinIdx{s}, subFreqLabels{s}, subBinLabels{s}] = selectDataForTraining(sourceDataFileName, rca_bins, rca_harms, rca_conds, rca_trials);
        train_data(:,s) = all_data(:,s);
    end

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

if compare_ch > nChannels
    error_msg = sprintf('\n The requested comparision channel (%d) exceeds the maximum number of channels in the dataset (%d), ignoring request...\n', compare_ch, nChannels);
    error(error_msg);
else
end

%% REORDER WITH BINS
if (length(unique(binIndices))) > 1 && (any(binIndices == 1)) 
    % only reorder data if more than one thing in the stack                                                         
    % and average is in the stack  
    % move averages to the end of the stack
    for b = 1:length(binIndices)
        if binIndices(b) == max(binIndices)
            new_idx(b) = b - (length(unique(binIndices))-1);
        else
            new_idx(b) = b + 1;
        end
    end
    new_idx = [new_idx, new_idx+length(new_idx)];
    all_data = cellfun(@(x) x(new_idx, :, :), all_data, 'uni', false);  
    % note: not necessary to reorder train_data, as order of components 
    % should not influence the output
else
end

% always reorder bin indices, labels and rca_bins
binIndices = binIndices(new_idx(1:length(binIndices)));
newIndices(binIndices == 1) = length(binLabels); % binLabels has the full set of bins
newIndices(binIndices > 1) = binIndices(binIndices > 1) - 1;
binIndices = reshape(newIndices, [], 1); 
binLabels = [binLabels(2:end); binLabels(1)];
% make sure rca bins = 0 refer to the now last index ("ave")
rca_bins(rca_bins == 0) = length(binLabels);
clear new_idx;


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
if isempty(rca_harms)
    rca_harms = 'all';
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
    if ~return_all
        % if not returning all, also restrict all data
        all_data = all_data(:, rca_subs);
    else
    end
else
    rca_subs = 'all';
end

%% run RCA
fprintf('Running RCA...\n');
warning('off','all')
[~, W, A, Rxx, Ryy, Rxy, dGen] = rcaRun(train_data, n_reg, n_comp); 
COV.Rxx = Rxx; 
COV.Ryy = Ryy;
COV.Rxy = Rxy;
COV.sortedGeneralizedEigenValues = dGen;
% rca was trained on 'train_data', 
% but now project all the data through those components
rcaData = rcaProject(all_data, W);
lo_noise = rcaProject(cellNoiseData1, W); 
hi_noise = rcaProject(cellNoiseData2, W);

%% CREATE A COMPONENT OF JUST ONE CHANNEL, FOR PERFORMANCE EVALUATION
W_comp = zeros(nChannels,1); W_comp(compare_ch) = 1; 
rcaData = cellfun(@(x, y) cat(2, x, y), rcaData, rcaProject(all_data, W_comp), 'uni', false);
lo_noise = cellfun(@(x, y) cat(2, x, y), lo_noise, rcaProject(cellNoiseData1, W_comp), 'uni', false);
hi_noise = cellfun(@(x, y) cat(2, x, y), hi_noise, rcaProject(cellNoiseData2, W_comp), 'uni', false);
compLabels = [arrayfun(@(x) sprintf('rc%02d',x), 1:n_comp, 'uni', false), sprintf('channel%d', compare_ch)];
warning('on','all')

%% package some of the necessary variables for plotting/grouping data later on
rcaSettings.rcaBins = rca_bins;
rcaSettings.rcaFreqs = rca_harms;
rcaSettings.rcaConds = rca_conds;
rcaSettings.rcaTrials = rca_trials;
rcaSettings.rcaSubs = rca_subs;
rcaSettings.n_reg = n_reg;
rcaSettings.n_comp = n_comp;
rcaSettings.dataFolders = path_names;
rcaSettings.freqLabels = freqLabels; 
rcaSettings.binLabels = binLabels;
rcaSettings.compLabels = compLabels;
rcaSettings.data_type = data_type;
runDate = datestr(clock,26);
runDate(strfind(runDate,'/')) ='';
rcaSettings.runDate = runDate; % on what date was this RCA run?

%% generate final output struct array, organized by frequency in input data

%% THIS SECTION NEEDS TO BE MODERATED
% put shared values across frequencies into output struct
% this includes empty fields "mean", "subjects", and "stats"
% which will be used for aggregating data later (see aggregateData fxn)
rca_struct = struct('W', W, 'A', A, 'COV', COV, ...
    'mean', struct([]), 'subjects', struct([]), 'stats', struct([]));

all_freqs = unique(freqIndices);
% make an array of structs
rca_struct = repmat(rca_struct, length(all_freqs), 1);
% and populate it
for f = 1:length(all_freqs)
    fIdx = repmat(freqIndices==all_freqs(f),2,1); % repmat because the first half is real, second half is imag with same ordering
    rca_struct(f).rca_data = cellfun(@(x) x(fIdx,:,:),rcaData,'uni',false);
    rca_struct(f).noiseLower = cellfun(@(x) x(fIdx,:,:), lo_noise, 'uni',false);
    rca_struct(f).noiseHigher = cellfun(@(x) x(fIdx,:,:), hi_noise, 'uni',false);
    rca_struct(f).input_data = cellfun(@(x) x(fIdx,:,:),all_data,'uni',false);
    % output indices that match the data
    rcaSettings.freqIndices = freqIndices(freqIndices==all_freqs(f));
    rcaSettings.binIndices = binIndices(freqIndices==all_freqs(f));
    rca_struct(f).settings = orderfields(rcaSettings);
end