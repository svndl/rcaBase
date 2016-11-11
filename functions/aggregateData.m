function [avgData,realSubjs,imagSubjs] = aggregateData(rcaData,rcaSettings,keepConditions,ampErrorType,trialError)
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
if (nargin<5 || isempty(trialError)), trialError=false; else end


nConditions = size(rcaData,1);
nSubjects = size(rcaData,2);
nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);
nCompFromInputData = max(max(cellfun(@(x) size(x,2),rcaData)));
nTrials = max(max(cellfun(@(x) size(x,3),rcaData)));

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
            if trialError
                % grab all subjects' data, without averaging over trials: 
                tempReal = [];
                tempImag = [];
                for s=1:nSubjects
                    tempReal = cat(3,tempReal,rcaDataReal{condNum,s});
                    tempImag = cat(3,tempImag,rcaDataImag{condNum,s});
                end
                muRcaDataRealAllSubj(:,:,:,condNum) = nan(nBins*nFreqs,nCompFromInputData,nSubjects*nTrials);
                muRcaDataImagAllSubj(:,:,:,condNum) = nan(nBins*nFreqs,nCompFromInputData,nSubjects*nTrials);
                muRcaDataRealAllSubj(:,:,1:size(tempReal,3),condNum) = tempReal; % 1:size(tempReal,3) - make it possible to have subjects with fewer trials
                muRcaDataImagAllSubj(:,:,1:size(tempImag,3),condNum) = tempImag;
            else
                % grab all subjects' data, averaging over trials: 
                muRcaDataRealAllSubj(:,:,:,condNum) = nan(nBins*nFreqs,nCompFromInputData,nSubjects);
                muRcaDataImagAllSubj(:,:,:,condNum) = nan(nBins*nFreqs,nCompFromInputData,nSubjects);
                for s=1:nSubjects
                    muRcaDataRealAllSubj(:,:,s,condNum) = nanmean(rcaDataReal{condNum,s},3);
                    muRcaDataImagAllSubj(:,:,s,condNum) = nanmean(rcaDataImag{condNum,s},3);
                end
            end
            
            for rc=1:nCompFromInputData
                for f=1:nFreqs
                    for b=1:nBins
                        firstDimInd = rcaSettings.freqIndices==rcaSettings.freqsToUse(f) & rcaSettings.binIndices==rcaSettings.binsToUse(b);
                        realSubjs(b,:,f,rc,condNum) = squeeze(muRcaDataRealAllSubj(firstDimInd,rc,:,condNum));
                        imagSubjs(b,:,f,rc,condNum) = squeeze(muRcaDataImagAllSubj(firstDimInd,rc,:,condNum));
                        xyData = [realSubjs(b,:,f,rc,condNum)' imagSubjs(b,:,f,rc,condNum)'];
                        if size(xyData,1)<2
                            keyboard;
                        end
                        nanVals = sum(isnan(xyData),2)>0;                        
                        ampErrBins(b,f,rc,condNum,:) = fitErrorEllipse(xyData(~nanVals,:),ampErrorType);
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
    
    %     if calcErrors
    %         withinSubj = false;
    %         ampErrBins = nan(nBins,nFreqs,nCompFromInputData,2);
    %         % grab all subjects' data, averaging over trials:
    %         muRcaDataRealAllSubj = nanmean(cat(4,rcaDataReal{:}),3);
    %         muRcaDataImagAllSubj = nanmean(cat(4,rcaDataImag{:}),3);
    %         for rc=1:nCompFromInputData
    %             for f=1:nFreqs
    %                 for b=1:nBins
    %                     firstDimInd = rcaSettings.freqIndices==rcaSettings.freqsToUse(f) & rcaSettings.binIndices==rcaSettings.binsToUse(b);
    %                     realSubjs = squeeze(muRcaDataRealAllSubj(firstDimInd,rc,:));
    %                     imagSubjs = squeeze(muRcaDataImagAllSubj(firstDimInd,rc,:));
    %                     xyData = [realSubjs imagSubjs];
    %                     [~,ampErrBins(b,f,rc,:)] = getErrorEllipse(xyData,withinSubj,ampErrorType);
    %                 end
    %             end
    %         end
    %     end
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
