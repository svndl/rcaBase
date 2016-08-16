function [subjData] = getSubjectData(rcaData,rcaSettings,compNum,condNum,freqNum)
% [subjData] = getSubjectData(rcaData,rcaSettings,compNum,condNum,freqNum)
%
% return a nBins x 2 x nSubjects matrix of the real and imaginary
% components of the requested reliable component, compNum, for condition
% condNum at frequency freqNum.
%
% This matrix is useful for constructing the nBins x 6 x nSubjects matrix
% used by the scoring algorithm, as well as for analyzing effects across
% subjects.

if nargin<2, error('You must provide the rcaSettings struct created when your rca data were created.'); end

nConditions = size(rcaData,1);
nSubjects = size(rcaData,2);
nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);
nCompFromInputData = size(rcaData{1,1},2);

% convert to real/imaginary
[rcaDataReal,rcaDataImag] = getRealImag(rcaData);

subjData = nan(nBins,2,nSubjects);

rc = compNum;
f = freqNum;
for s = 1:nSubjects
    subjData(:,1,s) = nanmean(rcaDataReal{condNum,s}(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,:),3);
    subjData(:,2,s) = nanmean(rcaDataImag{condNum,s}(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,:),3);
end

