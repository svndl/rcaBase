function [avgData] = aggregateData(rcaData,rcaSettings,keepConditions,ampErrorType)
% [avgData] = aggregateBins(rcaData,rcaSettings,[keepConditions],[ampErrorType])
%
% rcaData: created during call to rcaSweep
% rcaSettings: created during call to rcaSweep
%
% avgData is a struct containing subject and trial averaged data that
% contains the following fields:
%
%   realBins: bin-by-harmonic-by-component array of REAL RC coefficients
%   imagBins: bin-by-harmonic-by-component array of IMAGINARY RC coefficients
%   ampBins: bin-by-harmonic-by-component array of RC amplitudes
%   phaseBins: bin-by-harmonic-by-component array of RC phases
%
%
% if keepConditions is true (default = false), condition number adds a 4th
%   dimension to each array.
% if ampErrorType is specified, then an additional field is returned:
%
%   ampErrBins: bin-by-harmonic-by-component array of RC amplitude errors
%       as estimated by specified ampErrorType: 'SEM' '95CI' or a string specifying
%       a different percentage CI formated following: '%.1fCI'. Default is
%       'SEM' (uses function getErrorEllipse.m).
%%
if nargin<2, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if nargin<3, keepConditions=false; end
if (nargin<4 || isempty(ampErrorType)), calcErrors=false; else calcErrors = true; end

nConditions = size(rcaData,1);
nSubjects = size(rcaData,2);
nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);
nCompFromInputData = size(rcaData{1,1},2);

% convert to real/imaginary
[rcaDataReal,rcaDataImag] = getRealImag(rcaData);

if keepConditions 
    realBins = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    imagBins = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    muRcaDataReal = nan(nBins*nFreqs,nCompFromInputData,nConditions);
    muRcaDataImag = nan(nBins*nFreqs,nCompFromInputData,nConditions);       
    if calcErrors        
        ampErrBins = nan(nBins,nFreqs,nCompFromInputData,nConditions,2);
    end
    
    for condNum = 1:nConditions        
        % average over trials and subjects for each condition
        muRcaDataReal(:,:,condNum) = nanmean(cat(3,rcaDataReal{condNum,:}),3); % weights each trial equally
        muRcaDataImag(:,:,condNum) = nanmean(cat(3,rcaDataImag{condNum,:}),3); % weights each trial equally
        
        for rc=1:nCompFromInputData
            for f=1:nFreqs
                realBins(:,f,rc,condNum) = muRcaDataReal(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,condNum);
                imagBins(:,f,rc,condNum) = muRcaDataImag(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,condNum);
            end
        end     
        
        if calcErrors
            % grab all subjects' data, averaging over trials: 
            muRcaDataRealAllSubj = nan(nBins*nFreqs,nCompFromInputData,nConditions);
            muRcaDataImagAllSubj = nan(nBins*nFreqs,nCompFromInputData,nConditions);
            for s=1:nSubjects
                muRcaDataRealAllSubj(:,:,s) = nanmean(rcaDataReal{condNum,s},3);
                muRcaDataImagAllSubj(:,:,s) = nanmean(rcaDataImag{condNum,s},3);
            end
            
            for rc=1:nCompFromInputData
                for f=1:nFreqs
                    for b=1:nBins
                        firstDimInd = rcaSettings.freqIndices==rcaSettings.freqsToUse(f) & rcaSettings.binIndices==rcaSettings.binsToUse(b);
                        realSubjs = squeeze(muRcaDataRealAllSubj(firstDimInd,rc,:));
                        imagSubjs = squeeze(muRcaDataImagAllSubj(firstDimInd,rc,:));
                        xyData = [realSubjs imagSubjs];
                        if size(xyData,1)<2
                            keyboard;
                        end
                        ampErrBins(b,f,rc,condNum,:) = fitErrorEllipse(xyData,ampErrorType);
                    end
                end
            end            
        end
    end
    
else    
    % average over trials, subjects AND conditions
    muRcaDataReal = nanmean(cat(3,rcaDataReal{:}),3);
    muRcaDataImag = nanmean(cat(3,rcaDataImag{:}),3);
    
    realBins = zeros(nBins,nFreqs,nCompFromInputData);
    imagBins = zeros(nBins,nFreqs,nCompFromInputData);
    for rc = 1:nCompFromInputData
        for f = 1:nFreqs
            realBins(:,f,rc)=muRcaDataReal(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc);
            imagBins(:,f,rc)=muRcaDataImag(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc);            
        end
    end
    
    if calcErrors
        withinSubj = false;
        ampErrBins = nan(nBins,nFreqs,nCompFromInputData,2);
        % grab all subjects' data, averaging over trials:
        muRcaDataRealAllSubj = nanmean(cat(4,rcaDataReal{:}),3);
        muRcaDataImagAllSubj = nanmean(cat(4,rcaDataImag{:}),3);
        for rc=1:nCompFromInputData
            for f=1:nFreqs
                for b=1:nBins
                    firstDimInd = rcaSettings.freqIndices==rcaSettings.freqsToUse(f) & rcaSettings.binIndices==rcaSettings.binsToUse(b);
                    realSubjs = squeeze(muRcaDataRealAllSubj(firstDimInd,rc,:));
                    imagSubjs = squeeze(muRcaDataImagAllSubj(firstDimInd,rc,:));
                    xyData = [realSubjs imagSubjs];
                    [~,ampErrBins(b,f,rc,:)] = getErrorEllipse(xyData,withinSubj,ampErrorType);
                end
            end
        end
    end
end

ampBins = sqrt(realBins.^2+imagBins.^2);
phaseBins = atan(imagBins./realBins);

avgData.realBins = realBins;
avgData.imagBins = imagBins;
avgData.ampBins = ampBins;
avgData.phaseBins = phaseBins;

if calcErrors
    avgData.ampErrBins = ampErrBins;
    avgData.ampErrType = ampErrorType;
end
