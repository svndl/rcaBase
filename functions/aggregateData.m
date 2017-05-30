function [avgData,realSubjs,imagSubjs] = aggregateData(rcaData,rcaSettings,keepConditions,ampErrorType,trialError,nrFit)
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
    
    %   nrFit: freq x rc x condition logical, indicating whether or not to do Naka-Rushton fitting

    %% 

    if nargin<2, error('You must provide the rcaSettings struct created when your rca data were created.'); end
    if nargin<3, keepConditions=false; end
    if (nargin<4 || isempty(ampErrorType)), calcErrors=false; else calcErrors = true; end
    if (nargin<5 || isempty(trialError)), trialError=false; else end
    if (nargin<6 ), nrFit=[]; else end

    nSubjects = size(rcaData,2);

    if ~keepConditions
        % concatenate conditions, but keep everything else seperate
        rcaData = arrayfun(@(x) cat(3,rcaData{:,x}),1:nSubjects,'uni',false);
    else
    end

    nConditions = size(rcaData,1);
    nFreqs = length(rcaSettings.freqsToUse);
    nBins = length(rcaSettings.binsToUse);
    nCompFromInputData = max(max(cellfun(@(x) size(x,2),rcaData)));
    nTrials = max(max(cellfun(@(x) size(x,3),rcaData)));

    % convert to real/imaginary
    [rcaDataReal,rcaDataImag] = getRealImag(rcaData);

    if calcErrors        
        ampErrBins = nan(nBins,nFreqs,nCompFromInputData,nConditions,2);
        tPval = nan(nBins,nFreqs,nCompFromInputData,nConditions);
        tSqrd = nan(nBins,nFreqs,nCompFromInputData,nConditions);
        tSig = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    end

    realBins = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    imagBins = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    muRcaDataReal = nan(nBins*nFreqs,nCompFromInputData,nConditions);
    muRcaDataImag = nan(nBins*nFreqs,nCompFromInputData,nConditions);
    
    NR_pOpt = nan(4,nFreqs,nCompFromInputData,nConditions);
    NR_R2 = nan(1,nFreqs,nCompFromInputData,nConditions);
    NR_Range = nan(1,nFreqs,nCompFromInputData,nConditions);
    NR_hModel = [];
    NR_pOpt_JKSE = nan(4,nFreqs,nCompFromInputData,nConditions);
    
    % assign default value to nrFit
    if isempty(nrFit)
        nrFit = false(nFreqs,nCompFromInputData,nConditions);
    else
    end
    
    for condNum = 1:nConditions
        % average over trials and subjects for each condition
        muRcaDataReal(:,:,condNum) = nanmean(cat(3,rcaDataReal{condNum,:}),3); % weights each trial equally
        muRcaDataImag(:,:,condNum) = nanmean(cat(3,rcaDataImag{condNum,:}),3); % weights each trial equally
        % make sure only selected components are included
        binLevels = cell2mat(cellfun(@(x) str2num(x), rcaSettings.binLevels{condNum},'uni',false));
        for rc=1:nCompFromInputData
            for f=1:nFreqs
                realBins(:,f,rc,condNum) = muRcaDataReal(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,condNum);
                imagBins(:,f,rc,condNum) = muRcaDataImag(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,condNum);
                 % fit Naka Rushton on means
                if nrFit(f,rc,condNum)
                    fitData = sqrt(realBins(:,f,rc,condNum).^2+imagBins(:,f,rc,condNum).^2);
                    [ NR_pOpt(:,f,rc,condNum), NR_R2(:,f,rc,condNum), NR_Range(:,f,rc,condNum), NR_hModel ] = FitNakaRushton(binLevels, fitData);
                else
                end
            end
        end
       
    end


    if calcErrors  
        for condNum = 1:nConditions
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
                        % make sure only selected components are included
                        firstDimInd = rcaSettings.freqIndices==rcaSettings.freqsToUse(f) & rcaSettings.binIndices==rcaSettings.binsToUse(b);
                        realSubjs(b,:,f,rc,condNum) = squeeze(muRcaDataRealAllSubj(firstDimInd,rc,:,condNum));
                        imagSubjs(b,:,f,rc,condNum) = squeeze(muRcaDataImagAllSubj(firstDimInd,rc,:,condNum));
                        xyData = [realSubjs(b,:,f,rc,condNum)' imagSubjs(b,:,f,rc,condNum)'];
                        if size(xyData,1)<2
                            keyboard;
                        end
                        nanVals = sum(isnan(xyData),2)>0;                        
                        ampErrBins(b,f,rc,condNum,:) = fitErrorEllipse(xyData(~nanVals,:),ampErrorType);
                        % compute t2-statistic against zero
                        tStruct = tSquaredFourierCoefs(xyData(~nanVals,:));
                        tPval(b,f,rc,condNum) = tStruct.pVal;
                        tSqrd(b,f,rc,condNum) = tStruct.tSqrd;
                        tSig(b,f,rc,condNum) = tStruct.H;
                    end
                    if nrFit(f,rc,condNum)
                        for s = 1:size(realSubjs,2) % number of subjects, or subject x trials
                            sIdx = true(size(realSubjs,2),1);
                            sIdx(s) = false;
                            fitData(:,s) = sqrt( nanmean(realSubjs(:,sIdx,f,rc,condNum),2).^2 + nanmean(imagSubjs(:,sIdx,f,rc,condNum),2).^2 );
                        end
                        NR_JK_pOpt = FitNakaRushton(binLevels, fitData);
                        NR_pOpt_JKSE(:,f,rc,condNum) = getParSE( NR_JK_pOpt );
                        clear NR_JK_pOpt;
                    else
                    end
                end
            end
        end
    else
    end

    ampBins = sqrt(realBins.^2+imagBins.^2);
    phaseBins = atan(imagBins./realBins);

    avgData.realBins = realBins;
    avgData.imagBins = imagBins;
    avgData.ampBins = ampBins;
    avgData.phaseBins = phaseBins;

    % Naka-Rushton output
    avgData.NakaRushton.pOpt = NR_pOpt;
    avgData.NakaRushton.R2 = NR_R2;
    avgData.NakaRushton.Range = NR_Range;
    avgData.NakaRushton.hModel = NR_hModel;
    
    if calcErrors
        avgData.ampErrBins = ampErrBins;
        avgData.ampErrType = ampErrorType;
        avgData.tSqrdSig = tSig;
        avgData.tSqrdP = tPval;
        avgData.tSqrdVal = tSqrd;
        avgData.NakaRushton.JKSE = NR_pOpt_JKSE;
    end
end

function pSE = getParSE( pJK )
        % function from Spero for computing jack-knifed SEs for NR parameters
		nDim = ndims( pJK );
		nSubj = size( pJK, nDim );
		pSE = sqrt( ( nSubj - 1 ) / nSubj * sum( bsxfun( @minus, pJK, mean(pJK,nDim) ).^2, nDim ) );
% 		p(:) = nSubj * p - ( nSubj - 1 ) * mean( pJK, nDim );			% unbiased estimate of mean.  ???making things go negative???, bar graphs should match plot anyhow
end
