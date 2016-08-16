function [rcaData,W,A,covData,noiseData,comparisonData,comparisonNoiseData,rcaSettings]=rcaSweep(pathnames,binsToUse,freqsToUse,condsToUse,nReg,nComp,dataType,chanToCompare,show,rcPlotStyle)
% perform RCA on sweep SSVEP data exported to RLS or DFT format
%
% [rcaData,W,A,noiseData,compareData,compareNoiseData,freqIndices,binIndices]=RCASWEEP(PATHNAMES,[BINSTOUSE],[FREQSTOUSE],[CONDSTOUSE],[NREG],[NCOMP],[DATATYPE],[COMPARECHAN],[SHOW],[RCPLOTSTYLE])
%
% INPUTS:
% pathnames (required): cell vector of string directory names housins DFT_c00x.txt or RLS_c00x.txt exports
%   for example,  pathnames={'/Volumes/Denali_4D2/rca/s001/','/Volumes/Denali_4D2/rca/s002/','/Volumes/Denali_4D2/rca/s003/'}
% binsToUse: vector of bin indices to include in RCA (defaults to bin 0 or average across bins)
% freqsToUse: vector of frequency indices to include in RCA (defaults to 1 or first harmonic ?)
% condsToUse: vector of conditions to use
% nReg: RCA regularization parameter (defaults to 9)
% nComp: number of RCs to retain (defaults to 3)
% dataType: can be 'DFT' or 'RLS'
% compareChan: comparison channel index between 1 and the total number
%   of channels in the specified dataset
% show: 1 to see a figure of sweep amplitudes for each harmonic and component (defaults to 1), 0 to not display
% rcPlotStyle: see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
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

if nargin<10 || isempty(rcPlotStyle), rcPlotStyle = []; end
if nargin<9 || isempty(show), show=1; end
if nargin<8 || isempty(chanToCompare), 
    computeComparison=false; 
    chanToCompare=NaN;
elseif ~isempty(chanToCompare)
    computeComparison=true; 
end
if nargin<7 || isempty(dataType), dataType = 'RLS'; end
if nargin<6 || isempty(nComp), nComp=3; end
if nargin<5 || isempty(nReg), nReg=9; end
if nargin<4 || isempty(condsToUse), condsToUse=1; end
if nargin<3 || isempty(freqsToUse), freqsToUse=1; end
if nargin<2 || isempty(binsToUse), binsToUse=0; end
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
for s=1:nSubjects
    sourceDataFileName = sprintf('%s/sourceData_%s.mat',pathnames{s},dataType);
    if isempty(dir(sourceDataFileName))
        createSourceDataMat(pathnames{s});       
    end        
    [signalData,indF,indB,noise1,noise2,freqLabels,binLevels] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse);
    freqIndices{s}=indF;
    binIndices{s}=indB;
    sensorData(:,s)=signalData;
    cellNoiseData1(:,s)=noise1;
    cellNoiseData2(:,s)=noise2;
    fprintf('Done selecting data for subject %d/%d: %s\n',s,nSubjects,pathnames{s});
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

%% run RCA
fprintf('Running RCA...\n');
[rcaData,W,A,Rxx,Ryy,Rxy,dGen,plotSettings]=rcaRun(sensorData,nReg,nComp,[],[],show,rcPlotStyle); 
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

%% package some of the necessary variables for plotting/grouping data later on
rcaSettings.dataFolders = pathnames;
rcaSettings.freqIndices = freqIndices;
rcaSettings.binIndices = binIndices;
rcaSettings.binsToUse = binsToUse;
rcaSettings.freqsToUse = freqsToUse;
rcaSettings.condsToUse = condsToUse;
rcaSettings.nReg = nReg;
rcaSettings.nComp = nComp;
rcaSettings.chanToCompare = chanToCompare;
rcaSettings.freqLabels = freqLabels; 
rcaSettings.binLevels = binLevels;
rcaSettings.dataType = dataType;
rcaSettings.RCplottingInfo = plotSettings;