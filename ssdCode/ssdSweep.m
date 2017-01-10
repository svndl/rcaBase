function ssdStrct=ssdSweep(pathnames,binsToUse,freqsToUse,condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,show,rcPlotStyle,doReg)
% perform SSD on sweep SSVEP data exported to RLS or DFT format
% different to rca, training is done in the complex domain
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
% show: 1 to see a figure of sweep amplitudes for each harmonic and component (defaults to 1), 0 to not display
% rcPlotStyle: see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
%
% OUTPUTS:
% ssdStruct: returned by ssdRun
% W: linear transformation matrix to go from sensor space to RC-space
% A: linear transformation matrix to go from RC-space space to sensor space
% covData: the covariance data used to learn the RCs
% noiseData: struct containing fields for the rca space representation of frequency side bands
% compareData: rcaData for the comparison channel only (if supplied)
% compareNoiseData: noiseData for the comparison channel only (if supplied)
% ssdSettings: struct containing settings used to run RCA & indices
%   necessary for unpacking rcaData later on
%
% ssdData is a (NumberConditions x NumberSubjects) cell array. Each cell
% contains a (length(freqsToUse)*length(binsToUse)*2 x nComp x NumberTrialsPerSubject)
% matrix containing the real and imaginary components of each RC.
%
% file adapted from rcaSweep by SB, 12/2016

if nargin<11 || isempty(rcPlotStyle), rcPlotStyle = []; end
if nargin<10 || isempty(show), show=1; end
if nargin<9 || isempty(chanToCompare), 
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
freqIndices=cell(nSubjects,1);
binIndices=cell(nSubjects,1);
fprintf('Reading in sensor data from provided path names...\n')
s = 1;
while (s <= nSubjects)
    sourceDataFileName = sprintf('%s/sourceData_%s.mat',pathnames{s},dataType);
    if isempty(dir(sourceDataFileName))
        createSourceDataMat(pathnames{s});
    end
    [signalData,freqIndices{s},binIndices{s},noise1,noise2,freqLabels,binLevels] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse);
    sensorData(:,s)=signalData;
    cellNoiseData1(:,s)=noise1;
    cellNoiseData2(:,s)=noise2;
    fprintf('Done selecting data for subject %d/%d: %s\n',s,nSubjects,pathnames{s});
    s = s + 1;
end
nChannels = size(sensorData{1,1},2);
if computeComparison
    if chanToCompare>nChannels
        fprintf('The requested comparision channel (%d) exceeds the maximum number of channels in the dataset (%d), ignoring request...\n',chanToCompare,nChannels)
        computeComparison = false;
    end
end

%% check freq and bin indices for consistency across subjects
for s=1:nSubjects
    if sum(abs(freqIndices{s}-freqIndices{1}))~=0 && sum(abs(binIndices{s}-binIndices{1}))~=0
        error('Frequency and bin indices vary across subjects: check consistency of DFT/RLS exports\n.');
    end
end
% if they're all the same, make into a regular array
freqIndices=freqIndices{1};
binIndices=binIndices{1};
%% run SSD

fprintf('Running SSD...\n');
[ssdData,W,A,C_s,C_n,dGen,plotSettings]=ssdRun_rcaFormat(sensorData,cellNoiseData1,cellNoiseData2,nReg,nComp,[],[],show,rcPlotStyle); 

covData.C_s = C_s;
covData.C_n = C_n;
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

%% package some of the necessary variables for plotting/grouping data later on
rcaSettings.dataFolders = pathnames;
rcaSettings.freqIndices = freqIndices;
rcaSettings.binIndices = binIndices;
rcaSettings.binsToUse = binsToUse;
rcaSettings.freqsToUse = freqsToUse;
rcaSettings.condsToUse = condsToUse;
rcaSettings.trialsToUse = trialsToUse;
rcaSettings.nReg = nReg;
rcaSettings.nComp = nComp;
rcaSettings.chanToCompare = chanToCompare;
rcaSettings.freqLabels = freqLabels; 
rcaSettings.binLevels = binLevels;
rcaSettings.dataType = dataType;
rcaSettings.RCplottingInfo = plotSettings;

%% generate final output struct
ssdStrct.data = ssdData;
ssdStrct.W = W;
ssdStrct.A = A;
ssdStrct.covData = covData; 
ssdStrct.noiseData = noiseData;
ssdStrct.comparisonData = comparisonData;
ssdStrct.comparisonNoiseData = comparisonNoiseData;
ssdStrct.inputData = sensorData;
ssdStrct.settings = rcaSettings;