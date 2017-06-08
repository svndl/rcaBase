function [signalDataSel,indFSel,indBSel,noise1Sel,noise2Sel,freqLabelsSel,binLevelsSel,trialsSel] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse)
    %
    % [signalData,indF,indB,noise1,noise2,freqLabelsSel,binLevelsSel,trialsSel] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse)
    %
    % Return only the desired subset of data contained in sourceDataFileName
    % which contains only binsToUse, freqsToUse, condsToUse, trialsToUse
    load(sourceDataFileName); % should contain these variables: signalData,indF,indB,noise1,noise2,freqLabels,binLevels,chanIncluded

    if nargin < 5; trialsToUse = []; else end
    if nargin < 4 || isempty(condsToUse); 
        condsToUse = 1:size(signalData,1); 
    else
    end;
    if nargin < 3 
        freqsToUse = [];
    else
    end
    if nargin < 2
        binsToUse = [];
    else
    end

    nC = size(signalData,2);
    nR = length(condsToUse);
    signalDataSel = cell(nR,nC);
    noise1Sel = cell(nR,nC);
    noise2Sel = cell(nR,nC);
    binLevelsSel = cell(nR,1);
    trialsSel = cell(nR,1);
    
    for row = 1:nR
        if condsToUse(row) > size(signalData,1) || isempty(signalData{condsToUse(row)});
            signalDataSel(row,:) = {[]};
            noise1Sel(row,:) = {[]};
            noise2Sel(row,:) = {[]};
            binLevelsSel(row) = {[]};
            trialsSel(row) = {[]};
        else
            if isempty(binsToUse)
                binsToUse = unique(indB{condsToUse(row)});
                binsToUse = binsToUse(binsToUse>0); % use all bins except the average bin
            else
            end
            if isempty(freqsToUse)
                freqsToUse = unique(indF{condsToUse(row)});
            else
            end
            selRowIx = ismember(indB{condsToUse(row)},binsToUse) & ismember(indF{condsToUse(row)},freqsToUse);
            for col = 1:nC
                if isempty(trialsToUse)
                    trialsSel{row} = 1:size(signalData{condsToUse(row),col},3); % use all trials
                else
                    trialsSel{row} = trialsToUse;
                end
                missingIdx = find(~ismember(trialsSel{row},1:size(signalData{condsToUse(row),col},3)));
                if ~isempty(missingIdx)
                    error('Input trial indices "%s" is not among set of trials',num2str(trialsSel{row}(missingIdx),'%d,'));
                else
                end
                signalDataSel{row,col} = signalData{condsToUse(row),col}(repmat(selRowIx,[2,1]),:,trialsSel{condsToUse(row)}); % repmat because the first half is real, second half is imag with same ordering
                noise1Sel{row,col}     =     noise1{condsToUse(row),col}(repmat(selRowIx,[2,1]),:,trialsSel{condsToUse(row)});
                noise2Sel{row,col}     =     noise2{condsToUse(row),col}(repmat(selRowIx,[2,1]),:,trialsSel{condsToUse(row)});

                binLevelsSel{row} = binLevels{condsToUse(row)}(binsToUse+1); % add one, because bin level 0 = average
            end
            % check if indices are unequal
            nonEmpty = find(cell2mat(cellfun(@(x) ~isempty(x),indF(condsToUse),'uni',false)),1,'first');
            if any ( indF{nonEmpty} ~= indF{condsToUse(row)} )
                error('frequency indices are not matched across conditions');
            elseif any ( indB{nonEmpty} ~= indB{condsToUse(row)} )
                error('bin indices are not matched across conditions');
            else
            end
            if ~exist('indBSel','var') % bins and frequencies should not vary across conditions
                indBSel = indB{condsToUse(row)}(selRowIx);
                indFSel = indF{condsToUse(row)}(selRowIx);
            else
            end
        end
    end
    for k = 1:length(freqsToUse)
        freqLabelsSel{k,1} = freqLabels{freqsToUse(k)};
    end
    
    % now replace missing conditions with NaNs, dude!
    signalDataSel = replaceEmpty(signalDataSel);
    noise1Sel = replaceEmpty(noise1Sel);
    noise2Sel = replaceEmpty(noise2Sel);
    binLevelsSel = replaceEmpty(binLevelsSel);
    trialsSel = replaceEmpty(trialsSel);
end

function cellOut = replaceEmpty(cellIn)
    nonEmpty = cell2mat(cellfun(@(x) ~isempty(x),cellIn,'uni',false));
    repIdx = find(nonEmpty,1,'first');
    cellOut = cellIn;
    cellOut(~nonEmpty) = {nan(size(cellIn{repIdx}))};
end