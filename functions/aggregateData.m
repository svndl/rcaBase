function [avgData,muRcaDataRealAllSubj,muRcaDataImagAllSubj] = aggregateData(rcaStruct,keepConditions,ampErrorType,trialError,nrFit)
    % [avgData] = aggregateBins(rcaData,rcaSettings,[keepConditions],[ampErrorType])
    %
    % rcaStruct: created during call to rcaSweep
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
    % ampErrBins: bin-by-harmonic-by-component array of RC amplitude errors
    %      as estimated by specified ampErrorType: 'SEM' '95CI' or a string specifying
    %      a different percentage CI formated following: '%.1fCI'. Default is
    %      'SEM' (uses function getErrorEllipse.m).
    
    % nrFit: freq x rc x condition logical, indicating whether or not to do Naka-Rushton fitting

    %% 
    if nargin<2, keepConditions=false; end
    if (nargin<3 || isempty(ampErrorType)), ampErrorType = 'none'; else end
    if (nargin<4 || isempty(trialError)), trialError = false; else end
    if (nargin<5 ), nrFit=[]; else end
    
    rcaSettings = rcaStruct.settings;
    % add comparison data as last component
    if isfield(rcaStruct,'comparisonData')
        rcaData = cellfun(@(x,y) cat(2,x,y), rcaStruct.data,rcaStruct.comparisonData,'uni',false);
    else
        rcaData = rcaStruct.data;
    end
   
    nSubjects = size(rcaData,2);

    if ~keepConditions
        % concatenate conditions, but keep everything else seperate
        rcaData = arrayfun(@(x) cat(3,rcaData{:,x}),1:nSubjects,'uni',false);
    else
    end

    nConditions = size(rcaData,1);
    nFreqs = length(rcaSettings.freqsToUse);
    nBins = length(rcaSettings.binsToUse);
    nCompFromInputData = max(max(cellfun(@(x) size(x,2),rcaData))); % note, includes comparison data
    nTrials = max(max(cellfun(@(x) size(x,3),rcaData)));

    % do the noise
    % add comparison data, and compute the amplitudes, which is all you need
    if isfield(rcaStruct,'noiseData')
        ampNoiseBins = zeros(nBins,nFreqs,nCompFromInputData,nConditions);
        if trialError
            ampNoiseBinsSubjects = zeros(nBins,nFreqs,nCompFromInputData,nSubjects*nTrials,nConditions);
        else
            ampNoiseBinsSubjects = zeros(nBins,nFreqs,nCompFromInputData,nSubjects,nConditions);
        end
        for z = 1:2
            if z == 1
                % lower
                noiseStruct.data = rcaStruct.noiseData.lowerSideBand;
                noiseStruct.comparisonData = rcaStruct.comparisonNoiseData.lowerSideBand;
            else
                noiseStruct.data = rcaStruct.noiseData.higherSideBand;
                noiseStruct.comparisonData = rcaStruct.comparisonNoiseData.higherSideBand;
            end
            noiseStruct.settings = rcaSettings;
            noiseStruct.Out = aggregateData(noiseStruct,keepConditions,'none',trialError,[]); % do not compute NR or error
            ampNoiseBins = noiseStruct.Out.ampBins + ampNoiseBins;
            ampNoiseBinsSubjects = noiseStruct.Out.subjectAmp + ampNoiseBinsSubjects;
            clear noiseStruct;
        end
        ampNoiseBins = ampNoiseBins./2;
        ampNoiseBinsSubjects =  ampNoiseBinsSubjects./2;
    else
        ampNoiseBins = nan(nBins,nFreqs,nCompFromInputData,nConditions);
        ampNoiseBinsSubjects = nan(nBins,nFreqs,nCompFromInputData,nSubjects,nConditions);
    end
    
    % convert to real/imaginary
    [rcaDataReal,rcaDataImag] = getRealImag(rcaData);

    if ~strcmp(ampErrorType,'none')        
        ampErrBins = nan(nBins,nFreqs,nCompFromInputData,nConditions,2);
        tPval = nan(nBins,nFreqs,nCompFromInputData,nConditions);
        tSqrd = nan(nBins,nFreqs,nCompFromInputData,nConditions);
        tSig = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    end
    
    muRcaDataReal = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    muRcaDataImag = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    
    NR_Params = nan(5,nFreqs,nCompFromInputData,nConditions);
    NR_R2 = nan(1,nFreqs,nCompFromInputData,nConditions);
    NR_Range = nan(1,nFreqs,nCompFromInputData,nConditions);
    NR_hModel = [];
    NR_JK_SE = nan(5,nFreqs,nCompFromInputData,nConditions);
    NR_JK_Params = nan(5,nSubjects,nFreqs,nCompFromInputData,nConditions);
    
    % assign default value to nrFit
    if isempty(nrFit)
        nrFit = false(nFreqs,nCompFromInputData,nConditions);
    else
    end
    
    for condNum = 1:nConditions
        tempReal = [];
        tempImag = [];
        tempZ = [];
        for s=1:nSubjects
            if trialError
                nanSet = nan(nBins,nFreqs,nCompFromInputData,nSubjects*nTrials);
                % grab all subjects' data, without averaging over trials: 
                tempReal = cat(3,tempReal,rcaDataReal{condNum,s});
                tempImag = cat(3,tempImag,rcaDataImag{condNum,s});
            else
                nanSet = nan(nBins,nFreqs,nCompFromInputData,nSubjects);
                % grab all subjects' data, averaging over trials: 
                tempReal = cat(3,tempReal,nanmean(rcaDataReal{condNum,s},3));
                tempImag = cat(3,tempImag,nanmean(rcaDataImag{condNum,s},3));
            end
            tempZ = cat(3,tempZ, computeZsnr(rcaDataReal{condNum,s},rcaDataImag{condNum,s}) );
        end
        muRcaDataRealAllSubj(:,:,:,:,condNum) = nanSet;
        muRcaDataImagAllSubj(:,:,:,:,condNum) = nanSet;
        zRcaDataAllSubj(:,:,:,:,condNum) = nanSet;
        % split bins and frequencies, and make sure only selected components are included
        binLevels = cell2mat(cellfun(@(x) str2num(x), rcaSettings.binLevels{condNum},'uni',false));
        for rc = 1:nCompFromInputData
            for f = 1:nFreqs
                for b = 1:nBins
                    curIdx = rcaSettings.freqIndices{condNum}==rcaSettings.freqsToUse(f) & rcaSettings.binIndices{condNum}==rcaSettings.binsToUse(b);
                    muRcaDataRealAllSubj(b,f,rc,1:size(tempReal(curIdx,rc,:),3),condNum) = tempReal(curIdx,rc,:);
                    muRcaDataImagAllSubj(b,f,rc,1:size(tempImag(curIdx,rc,:),3),condNum) = tempImag(curIdx,rc,:);
                    zRcaDataAllSubj(b,f,rc,1:size(tempZ(curIdx,rc,:),3),condNum) = tempZ(curIdx,rc,:);
                end
            end
        end
        % average over trials and subjects for each condition
        muRcaDataReal(:,:,:,condNum) = nanmean(muRcaDataRealAllSubj(:,:,:,:,condNum),4); % weights each trial equally
        muRcaDataImag(:,:,:,condNum) = nanmean(muRcaDataImagAllSubj(:,:,:,:,condNum),4); % weights each trial equally
        zRcaData(:,:,:,condNum) = nanmean(zRcaDataAllSubj(:,:,:,:,condNum),4); % weights each trial equally
    end
    % fit Naka Rushton on means for all conditions
    if any(nrFit(:))
        fitData = sqrt(muRcaDataReal.^2+muRcaDataImag.^2);
        fitNoise = repmat(nanmean(ampNoiseBins,4),[1,1,1,nConditions]); % use same noise across conditions
        % [ NR_Params, NR_R2, NR_Range, NR_hModel ] = FitNakaRushton(binLevels, fitData, ampNoiseBins);
        % [ NR_Params, NR_R2, NR_hModel ] = FitNakaRushtonFixed(binLevels, fitData,ampNoiseBins);
        [ NR_Params, NR_R2, NR_hModel ] = FitNakaRushtonEq(binLevels, fitData);
        % [ NR_Params, NR_R2, NR_hModel ] = FitNakaRushtonNoise(binLevels, fitData,fitNoise);

        clear fitData;
    else
    end
    
    if ~strcmp(ampErrorType,'none')  
        for condNum = 1:nConditions
            for rc=1:nCompFromInputData
                for f=1:nFreqs
                    testReal(1:nBins,:) = squeeze(muRcaDataRealAllSubj(:,f,rc,:,condNum));
                    testImag(1:nBins,:) = squeeze(muRcaDataImagAllSubj(:,f,rc,:,condNum));
                    testNoise(1:nBins,:) = squeeze(nanmean(ampNoiseBinsSubjects(:,f,rc,:,:),5));
                    for b=1:nBins
                        xyData = [testReal(b,:)' testImag(b,:)'];
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
                        for s = 1:size(testReal,2) % number of subjects, or subject x trials
                            sIdx = true(size(testReal,2),1);
                            sIdx(s) = false;
                            jkData(:,s) = sqrt( nanmean(testReal(:,sIdx),2).^2 + nanmean(testImag(:,sIdx),2).^2 );
                            jkNoise(:,s) = nanmean(testNoise(:,sIdx),2);
                        end
                        %jkParams = FitNakaRushton(binLevels, jkData,jkNoise);
                        jkParams = FitNakaRushtonEq(binLevels, jkData);
                        %jkParams = FitNakaRushtonNoise(binLevels, jkData,jkNoise);
                        NR_JK_SE(:,f,rc,condNum) = jackKnifeErr( jkParams' );
                        NR_JK_Params(:,:,f,rc,condNum) = jkParams;
                        clear jkParams; clear jkData;
                    else
                    end
                end
            end
        end
    else
    end

    avgData.realBins = muRcaDataReal;
    avgData.imagBins = muRcaDataImag;
    avgData.ampBins = sqrt(muRcaDataReal.^2+muRcaDataImag.^2);
    avgData.ampNoiseBins = ampNoiseBins;
    avgData.phaseBins = atan(muRcaDataImag./muRcaDataReal);
    avgData.zSNR.mean = zRcaData;
    avgData.zSNR.subj = zRcaDataAllSubj;
    avgData.subjectAmp = sqrt(muRcaDataRealAllSubj.^2+muRcaDataImagAllSubj.^2);
    avgData.subjectAmpNoise = ampNoiseBinsSubjects;

    % Naka-Rushton output
    avgData.NakaRushton.Params = NR_Params;
    avgData.NakaRushton.R2 = NR_R2;
    %avgData.NakaRushton.Range = NR_Range;
    avgData.NakaRushton.hModel = NR_hModel;
    
    if ~strcmp(ampErrorType,'none') 
        avgData.ampErrBins = ampErrBins;
        avgData.ampErrType = ampErrorType;
        avgData.tSqrdSig = tSig;
        avgData.tSqrdP = tPval;
        avgData.tSqrdVal = tSqrd;
        avgData.NakaRushton.JackKnife.SE = NR_JK_SE;
        avgData.NakaRushton.JackKnife.Params = NR_JK_Params;
    end
end

function outZ = computeZsnr(realVals,imagVals)
    % move trial dimension to first
    realVals = permute(realVals,[3,2,1]);
    imagVals = permute(imagVals,[3,2,1]);
    nReal = size( realVals );
    nImag = size( imagVals );
    if any(nReal ~= nImag)
        error('real and imaginary values are not matched');
    else
    end
    nSets = prod( nReal(2:end) );
    outZ = nan( [1, nReal(2:end)]);
    for z = 1:nSets
        xyData = [realVals(:,z),imagVals(:,z)];
        nanVals = sum(isnan(xyData),2)>0;
        % use standard deviation, to compute the zSNR
        if size( xyData(~nanVals,:) ) > 2
            [ampErr,~,zSNR] = fitErrorEllipse(xyData(~nanVals,:),'1STD',false);
        else
            zSNR = NaN;
        end
        outZ(:,z) = zSNR;
    end
     % move trial dimension back to third
    outZ = permute(outZ,[3,2,1]);
end
 
% function pSE = getParSE( pJK )
%         % function from Spero for computing jack-knifed SEs for NR parameters
% 		nDim = ndims( pJK );
% 		nSubj = size( pJK, nDim );
% 		pSE = sqrt( ( nSubj - 1 ) / nSubj * sum( bsxfun( @minus, pJK, mean(pJK,nDim) ).^2, nDim ) );
% % 		p(:) = nSubj * p - ( nSubj - 1 ) * mean( pJK, nDim );			% unbiased estimate of mean.  ???making things go negative???, bar graphs should match plot anyhow
% end
