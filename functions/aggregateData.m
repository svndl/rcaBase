function rca_struct = aggregateData(rca_struct, keep_conditions, error_type, trial_wise, do_nr)
    % rca_struct = aggregateData(rca_struct, keep_conditions, error_type, trial_wise, do_nr)
    %
    % input: 
    %   rca_struct: created during call to rcaSweep
    %   keep_conditions: [true]/false
    %       if false, will average over conditions
    %   error_type: 'SEM'/'none'/'95CI'
    %       or a string specifyingn a different percentage 
    %       CI formated following: '%.1fCI'
    %   trial-wise: [true]/false
    %       if true, will compute trial-wise errors and stats
    %   do_nr: harmonic x rc x condition logical, 
    %       indicating when to do Naka-Rushton fitting
    %       (only makes sense to do for sweep data)
    %       default is not to do it anywhere
    %
    % output: function will return rca_struct, and add fields containing 
    %         aggregate statistics for all components (described below).
    %   
    %   Note that 
    %   (a) the comparison electrode data is automatically added as 
    %       as an additional component in the output.
    %   (b) if there is more than one bin in the input (sweep data) the 
    %       vector average across bins will be added as a final bin.
    %   (c) In the current iteration of rcaSweep, each rca_struct contains 
    %       only one harmonic, but harmonic is nonetheless preserved 
    %       as an output dimension
    %  
    %   "mean", with the following subfields:
    %       real_signal:        REAL coefficients, vector-averaged
    %       imag_signal:        IMAGINARY coefficients, vector-averaged
    %       amp_signal:         AMPLITUDES, vector-averaged
    %       phase_signal:       PHASES, vector-averaged
    %       z_snr:              Z_SNR
    %       amp_noise:          NOISE BAND AMPLITUDES, vector-averaged
    %   all subfields are bin-by-harmonic-by-component-by-condition arrays
    %
    %   "subjects", with the following subfields
    %       amp_signal:         AMPLITUDES
    %       proj_amp_signal:    PROJECTED AMPLITUDES 
    %       amp_noise:          NOISE BAND AMPLITUDES
    %       z_snr:              Z_SNR
    %   all subfields are bin-by-harmonic-by-component-by-subject-by-condition arrays
    %
    % "stats", with the following subfields
    %      amp_err_up: upper error bars
    %
    %      amp_err_lo: lower error bars
    %
    %      t_sig:      significance (true/false) for student's t-test against zero
    %                  run on projected amplitudes across subjects
    %           
    %      t_p:        p-values for student's t-test against zero
    %
    %      t_val:      t-values for student's t-test against zero
    %
    %      t2_sig:     significance (true/false) for Hotelling's T2 against zero
    %                  run on real and imag coefficients across subjects
    %           
    %      t2_p:       p-values for Hotelling's T2 against zero
    %
    %      t2_val:     t-values for Hotelling's T2 against zero
    %
    %      ... all of the above are
    %          bin-by-harmonic-by-component-by-condition arrays
    %
    %       err_type: string, what type of error bar was computed
    %
    %       naka_rushton: struct of naka-rushton outputs, pretty
    %              self-explanatory
    
    if nargin<2
        keep_conditions=true; 
    end
    if ( nargin<3 ) || isempty(error_type)
        error_type = 'SEM'; 
    else
    end
    if ( nargin<4 ) || isempty(trial_wise)
        trial_wise = false;
    else
    end
    if nargin<5
        do_nr=[]; 
    else
    end
    
    rcaSettings = rca_struct.settings;
    rcaData = rca_struct.rca_data;
   
    nSubjects = size(rcaData,2);

    if ~keep_conditions
        % concatenate conditions, but keep everything else seperate
        rcaData = arrayfun(@(x) cat(3,rcaData{:,x}),1:nSubjects,'uni',false);
    else
    end

    nConditions = size(rcaData,1);
    % get freqs and bins from the data
    dataFreqs = unique(rcaSettings.freqIndices);
    dataBins = unique(rcaSettings.binIndices);
    nFreqs = length(dataFreqs);
    nBins = length(dataBins);
    
    nTrials = max(max(cellfun(@(x) size(x,3),rcaData)));
    nCompFromInputData = max(max(cellfun(@(x) size(x,2),rcaData))); % note, includes comparison data
    nSampFromInputData = max(max(cellfun(@(x) size(x,1),rcaData)));
    if nFreqs*nBins*2 ~= nSampFromInputData
        error('number of frequencies and bins does not match the number of samples');
    else
    end
    
    % do the noise
    % add comparison data, and compute the amplitudes, which is all you need
    if isfield(rca_struct,'noiseLower')
        ampNoiseBins = zeros(nBins,nFreqs,nCompFromInputData,nConditions);
        if trial_wise
            ampNoiseBinsSubjects = zeros(nBins,nFreqs,nCompFromInputData,nSubjects*nTrials,nConditions);
        else
            ampNoiseBinsSubjects = zeros(nBins,nFreqs,nCompFromInputData,nSubjects,nConditions);
        end
        for z = 1:2
            if z == 1
                % lower
                noiseStruct.rca_data = rca_struct.noiseLower;
            else
                noiseStruct.rca_data = rca_struct.noiseHigher;
            end
            noiseStruct.settings = rcaSettings;
            noiseStruct.Out = aggregateData(noiseStruct, keep_conditions, 'none', []); % do not compute NR or error
            ampNoiseBins = noiseStruct.Out.mean.amp_signal + ampNoiseBins;
            ampNoiseBinsSubjects = noiseStruct.Out.subjects.amp_signal + ampNoiseBinsSubjects;
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

    if ~strcmp(error_type,'none')        
        ampErrBins = nan(nBins,nFreqs,nCompFromInputData,nConditions,2);
        tPval = nan(nBins,nFreqs,nCompFromInputData,nConditions);
        tSqrd = nan(nBins,nFreqs,nCompFromInputData,nConditions);
        tSig = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    end
    
    % nBins + 1, room for average
    muRcaDataReal = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    muRcaDataImag = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    
    NR_Params = nan(5,nFreqs,nCompFromInputData,nConditions);
    NR_R2 = nan(1,nFreqs,nCompFromInputData,nConditions);
    NR_Range = nan(1,nFreqs,nCompFromInputData,nConditions);
    NR_hModel = [];
    NR_JK_SE = nan(5,nFreqs,nCompFromInputData,nConditions);
    NR_JK_Params = nan(5,nSubjects,nFreqs,nCompFromInputData,nConditions);
    
    % assign default value to do_nr
    if isempty(do_nr)
        do_nr = false(nFreqs,nCompFromInputData,nConditions);
    else
    end
    
    for condNum = 1:nConditions
        tempReal = [];
        tempImag = [];
        tempZ = [];
        for s=1:nSubjects
            if trial_wise
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
        binLabels = cell2mat(cellfun(@(x) str2num(x), rcaSettings.binLabels,'uni',false));
        for rc = 1:nCompFromInputData
            for f = 1:nFreqs
                for b = 1:nBins
                    curIdx = find(rcaSettings.freqIndices==dataFreqs(f) & rcaSettings.binIndices==dataBins(b));
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
        
        % compute projected amplitudes
        realVector = permute(muRcaDataRealAllSubj(:,:,:,:,condNum),[4,1,2,3]);
        imagVector = permute(muRcaDataImagAllSubj(:,:,:,:,condNum),[4,1,2,3]);
        tempProject = vectorProjection(realVector,imagVector);
        projectAmpAllSubj(:,:,:,:,condNum) = permute(tempProject,[2,3,4,1]);
    end
    % fit Naka Rushton on means for all conditions
    if any(do_nr(:))
        % don't include average in fits
        fitData = sqrt(muRcaDataReal.^2+muRcaDataImag.^2);
        fitNoise = repmat(nanmean(ampNoiseBins,4),[1,1,1,nConditions]); % use same noise across conditions
        % [ NR_Params, NR_R2, NR_Range, NR_hModel ] = FitNakaRushton(binLabels, fitData, ampNoiseBins);
        % [ NR_Params, NR_R2, NR_hModel ] = FitNakaRushtonFixed(binLabels, fitData,ampNoiseBins);
        [ NR_Params, NR_R2, NR_hModel ] = FitNakaRushtonEq(binLabels, fitData(1:nBins,:,:,:));
        % [ NR_Params, NR_R2, NR_hModel ] = FitNakaRushtonNoise(binLabels, fitData,fitNoise);

        clear fitData;
    else
    end
    
    if ~strcmp(error_type,'none')  
        for condNum = 1:nConditions
            for rc=1:nCompFromInputData
                for f=1:nFreqs
                    % do include averages when computing error bars
                    testReal = reshape(muRcaDataRealAllSubj(:,f,rc,:,condNum), [], nSubjects);
                    testImag = reshape(muRcaDataImagAllSubj(:,f,rc,:,condNum), [], nSubjects);
                    for b=1:nBins
                        xyData = [testReal(b,:)' testImag(b,:)'];
                        if size(xyData,1)<2
                            keyboard;
                        end
                        nanVals = sum(isnan(xyData),2)>0;                        
                        ampErrBins(b,f,rc,condNum,:) = fitErrorEllipse(xyData(~nanVals,:),error_type);
                        % compute t2-statistic against zero
                        tStruct = tSquaredFourierCoefs(xyData(~nanVals,:));
                        tPval(b,f,rc,condNum) = tStruct.pVal;
                        tSqrd(b,f,rc,condNum) = tStruct.tSqrd;
                        tSig(b,f,rc,condNum) = tStruct.H;
                    end
                    % do not include averages when computing fitting errors
                    testReal = squeeze(muRcaDataRealAllSubj(1:(nBins-1),f,rc,:,condNum));
                    testImag = squeeze(muRcaDataImagAllSubj(1:(nBins-1),f,rc,:,condNum));
                    testNoise = nanmean(reshape(ampNoiseBinsSubjects(1:(nBins-1),f,rc,:,:),[],nSubjects, nConditions),3);
                    if do_nr(f,rc,condNum)
                        for s = 1:size(testReal,2) % number of subjects, or subject x trials
                            sIdx = true(size(testReal,2),1);
                            sIdx(s) = false;
                            jkData(:,s) = sqrt( nanmean(testReal(:,sIdx),2).^2 + nanmean(testImag(:,sIdx),2).^2 );
                            jkNoise(:,s) = nanmean(testNoise(:,sIdx),2);
                        end
                        %jkParams = FitNakaRushton(binLabels, jkData,jkNoise);
                        jkParams = FitNakaRushtonEq(binLabels, jkData);
                        %jkParams = FitNakaRushtonNoise(binLabels, jkData,jkNoise);
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
    
    % 
    mean_data.amp_signal = sqrt(muRcaDataReal.^2+muRcaDataImag.^2);
    mean_data.phase_signal = atan(muRcaDataImag./muRcaDataReal);
    mean_data.amp_noise = ampNoiseBins;
    mean_data.z_snr = zRcaData;
    rca_struct.mean = mean_data;
 
    sub_data.amp_signal = sqrt(muRcaDataRealAllSubj.^2+muRcaDataImagAllSubj.^2);
    sub_data.real_signal = muRcaDataRealAllSubj;
    sub_data.imag_signal = muRcaDataImagAllSubj;
    sub_data.amp_noise = ampNoiseBinsSubjects;
    sub_data.proj_amp_signal = projectAmpAllSubj;
    sub_data.z_snr = zRcaDataAllSubj;
    rca_struct.subjects = sub_data;
    
    % Naka-Rushton output
    stat_data.naka_rushton.Params = NR_Params;
    stat_data.naka_rushton.R2 = NR_R2;
    stat_data.naka_rushton.hModel = NR_hModel;
    
    if ~strcmp(error_type,'none') 
        % error bins
        stat_data.amp_lo_err = ampErrBins(:,:,:,:,1);
        stat_data.amp_up_err = ampErrBins(:,:,:,:,2);
        stat_data.err_type = error_type;
        % regular old t-test on projected amplitudes
        % one-tailed, would never be below zero
        t_params = {'alpha',0.05,'dim',4,'tail','right'};
        [stat_data.t_sig, stat_data.t_p, ci, sts] = ttest(rca_struct.subjects.proj_amp_signal,0,t_params{:});
        stat_data.t_sig = reshape(stat_data.t_sig, nBins, nFreqs, nCompFromInputData, []);
        stat_data.t_p = reshape(stat_data.t_p, nBins, nFreqs, nCompFromInputData, []);   
        stat_data.t_val = reshape(sts.tstat, nBins, nFreqs, nCompFromInputData, []);
        % Hotelling's T2
        stat_data.t2_sig = tSig;
        stat_data.t2_p = tPval;
        stat_data.t2_val = tSqrd;
        % Naka-Rushton
        stat_data.naka_rushton.JackKnife.SE = NR_JK_SE;
        stat_data.naka_rushton.JackKnife.Params = NR_JK_Params;
    end
    rca_struct.stats = stat_data;
    
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
