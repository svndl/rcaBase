function rca_struct = rcaSweep(pathnames,binsToUse,freqsToUse,condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,forceReload)
% perform RCA on sweep SSVEP data exported to RLS or DFT format
%
% [rcaData,W,A,noiseData,compareData,compareNoiseData,freqIndices,binIndices]=RCASWEEP(PATHNAMES,[BINSTOUSE],[FREQSTOUSE],[CONDSTOUSE],[TRIALSTOUSE],[NREG],[NCOMP],[DATATYPE],[COMPARECHAN],[SHOW],[RCPLOTSTYLE])
%
% INPUTS:
% pathnames (required): cell vector of string directory names housins DFT_c00x.txt or RLS_c00x.txt exports
%   for example,  pathnames={'/Volumes/Denali_4D2/rca/s001/','/Volumes/Denali_4D2/rca/s002/','/Volumes/Denali_4D2/rca/s003/'}
% binsToUse: vector of bin indices to include in RCA (defaults to bin 0 or average across bins)
% freqsToUse: vector of frequency indices to include in RCA (defaults to 1 or first harmonic ?)
% condsToUse: vector of conditions to use
% trailsToUse: vector of indices of trials to use
% nReg: RCA regularization parameter (defaults to 9)
% nComp: number of RCs to retain (defaults to 3)
% dataType: can be 'DFT' or 'RLS'
% compareChan: comparison channel index between 1 and the total number
%   of channels in the specified dataset
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
%   necessary for unpacking rcaData later on
%
% rcaData is a (NumberConditions x NumberSubjects) cell array. Each cell
% contains a (length(freqsToUse)*length(binsToUse)*2 x nComp x NumberTrialsPerSubject)
% matrix containing the real and imaginary components of each RC.
%
% Jacek P. Dmochowski, 2015, report bugs to dmochowski@gmail.com
% Edited by HEG 07/2015

if nargin<10 || isempty(forceReload), forceReload = false; end
if nargin<9 || isempty(chanToCompare)
    computeComparison=false; 
    chanToCompare=NaN;
elseif ~isempty(chanToCompare)
    computeComparison=true; 
end
if nargin<8 || isempty(dataType), dataType = 'RLS'; end
if nargin<7 || isempty(nComp), nComp=3; end
if nargin<6 || isempty(nReg), nReg=9; end
if nargin<5 || isempty(trialsToUse), trialsToUse = []; end
if nargin<4 || isempty(condsToUse), condsToUse=[]; end
if nargin<3 || isempty(freqsToUse), freqsToUse=[]; end
if nargin<2 || isempty(binsToUse), binsToUse=[]; end
if nargin<1, error('Must specify at least one input argument'); end

%% if pathanmes is a string, convert to cell
if ~iscell(pathnames)
    try
        pathnames=mat2cell(pathnames);
    catch
        error('unable to parse pathnames: check that it is a cell array of strings');
    end
end
if isempty(pathnames)
    error('Variable pathnames is empty.');
end
%% read in signal and noise data
nSubjects=numel(pathnames);
sensorData={};
cellNoiseData1={};
cellNoiseData2={};
subFreqIdx=cell(nSubjects,1);
subBinIdx=cell(nSubjects,1);
subFreqLabels=cell(nSubjects,1);
subBinLabels=cell(nSubjects,1);
fprintf('Reading in sensor data from provided path names...\n')
s = 1;
assigned = false(1,length(condsToUse));
while (s <= nSubjects)
    sourceDataFileName = sprintf('%s/sourceData_%s.mat',pathnames{s},dataType);
    if isempty(dir(sourceDataFileName)) || forceReload
        createSourceDataMat(pathnames{s});
    end
    [signalData,noise1,noise2,subFreqIdx{s},subBinIdx{s},subFreqLabels{s},subBinLabels{s}] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse);
    sensorData(:,s)=signalData;
    cellNoiseData1(:,s)=noise1;
    cellNoiseData2(:,s)=noise2;
    
    for c = 1:length(condsToUse)
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
    fprintf('Done selecting data for subject %d/%d: %s\n',s,nSubjects,pathnames{s});
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
nChannels = size(sensorData{1,1},2);
if computeComparison
    if chanToCompare>nChannels
        fprintf('The requested comparision channel (%d) exceeds the maximum number of channels in the dataset (%d), ignoring request...\n',chanToCompare,nChannels)
        computeComparison = false;
    end
end

%% run RCA
fprintf('Running RCA...\n');
warning('off','all')
[rcaData,W,A,Rxx,Ryy,Rxy,dGen]=rcaRun(sensorData,nReg,nComp); 
covData.Rxx = Rxx;
covData.Ryy = Ryy;
covData.Rxy = Rxy;
covData.sortedGeneralizedEigenValues = dGen;
noiseData.lowerSideBand=rcaProject(cellNoiseData1,W); 
noiseData.higherSideBand=rcaProject(cellNoiseData2,W);

%% create a "component" of just one channel for performance evaluation if requested
if computeComparison    
    wComparison=zeros(nChannels,1); wComparison(chanToCompare)=1; 
    comparisonData=rcaProject(sensorData,wComparison); 
    comparisonNoiseData.lowerSideBand =rcaProject(cellNoiseData1,wComparison); 
    comparisonNoiseData.higherSideBand =rcaProject(cellNoiseData2,wComparison); 
end
warning('on','all')

%% package some of the necessary variables for plotting/grouping data later on
rcaSettings.dataFolders = pathnames;
rcaSettings.rcaBins = binsToUse;
rcaSettings.rcaFreqs = freqsToUse;
rcaSettings.rcaConds = condsToUse;
rcaSettings.rcaTrials = trialsToUse;
rcaSettings.nReg = nReg;
rcaSettings.nComp = nComp;
rcaSettings.chanToCompare = chanToCompare;
rcaSettings.freqLabels = freqLabels; 
rcaSettings.binLabels = binLabels;
rcaSettings.dataType = dataType;
runDate = datestr(clock,26);
runDate(strfind(runDate,'/')) ='';
rcaSettings.runDate = runDate; % on what date was this RCA run?

%% generate final output struct array, organized by frequency in input data
for f = 1:length(freqsToUse)
    rca_struct(f).W = W;
    rca_struct(f).A = A;
    rca_struct(f).covData = covData;
    fIdx = repmat(freqIndices==freqsToUse(f),2,1); % repmat because the first half is real, second half is imag with same ordering
    rca_struct(f).data = cellfun(@(x) x(fIdx,:,:),rcaData,'uni',false);
    rca_struct(f).noiseData.lowerSideBand = cellfun(@(x) x(fIdx,:,:),noiseData.lowerSideBand,'uni',false);
    rca_struct(f).noiseData.higherSideBand = cellfun(@(x) x(fIdx,:,:),noiseData.higherSideBand,'uni',false);
    rca_struct(f).comparisonData = cellfun(@(x) x(fIdx,:,:),comparisonData,'uni',false);
    rca_struct(f).comparisonNoiseData.lowerSideBand = cellfun(@(x) x(fIdx,:,:),comparisonNoiseData.lowerSideBand,'uni',false);
    rca_struct(f).comparisonNoiseData.higherSideBand = cellfun(@(x) x(fIdx,:,:),comparisonNoiseData.higherSideBand,'uni',false);
    rca_struct(f).inputData = cellfun(@(x) x(fIdx,:,:),sensorData,'uni',false);
    % output indices that match the data
    rcaSettings.freqIndices = freqIndices(freqIndices==freqsToUse(f));
    rcaSettings.binIndices = binIndices(freqIndices==freqsToUse(f));
    % put into output struct
    rca_struct(f).settings = orderfields(rcaSettings);
    % fields for putting aggregate data later (see aggregateData fxn)
    rca_struct(f).mean = struct([]);
    rca_struct(f).subjects = struct([]); 
    rca_struct(f).stats = struct([]); 
end