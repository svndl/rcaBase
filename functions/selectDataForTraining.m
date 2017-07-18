function [signalDataSel,noise1Sel,noise2Sel,indFSel,indBSel,freqLabelsSel,binLabelsSel,trialsSel] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse)
    %
    % [signalData,noise1,noise2,indF,indB,,freqLabelsSel,binLevelsSel,trialsSel] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse)
    %
    % Return only the desired subset of data contained in sourceDataFileName
    % which contains only binsToUse, freqsToUse, condsToUse, trialsToUse
    load(sourceDataFileName); % should contain these variables: signalData,indF,indB,noise1,noise2,freqLabels,binLevels,chanIncluded

    if nargin < 5
        trialsToUse = []; 
    else
    end
    
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

    nConds = length(condsToUse);
    signalDataSel = cell(nConds,1);
    noise1Sel = cell(nConds,1);
    noise2Sel = cell(nConds,1);
    binLabelsSel = cell(nConds,1);
    freqLabelsSel = cell(nConds,1);
    trialsSel = cell(nConds,1);
    indBSel = cell(nConds,1);
    indFSel = cell(nConds,1);
    
    for c = 1:nConds
        if condsToUse(c) > size(signalData,1) || isempty(signalData{condsToUse(c)});
            signalDataSel{c} = [];
            noise1Sel{c} = [];
            noise2Sel{c} = [];
            indFSel{c} = [];
            indBSel{c} = [];
            freqLabelsSel{c} = [];
            binLabelsSel{c} = [];
            trialsSel{c} = [];
        else
            if isempty(binsToUse)
                % use all available bins
                % (note: canmot differ across conditions)
                binsToUse = unique(indB{condsToUse(c)});
                binsToUse = binsToUse(binsToUse>0); % use all bins except the average bin
            else
            end
            if isempty(freqsToUse)
                % use all available frequencies
                % (note: canmot differ across conditions)
                freqsToUse = unique(indF{condsToUse(c)});
            else
            end
            selRowIx = ismember(indB{condsToUse(c)},binsToUse) & ismember(indF{condsToUse(c)},freqsToUse);
            if isempty(trialsToUse)
                % use all available trials 
                % (note: can differ across conditions)
                curTrials = 1:size(signalData{condsToUse(c)},3); % use all trials
            else
                curTrials = trialsToUse;
            end
            missingIdx = find(~ismember(curTrials,1:size(signalData{condsToUse(c)},3)));
            if ~isempty(missingIdx)
                error('Input trial indices "%s" is not among set of trials',num2str(curTrials(missingIdx),'%d,'));
            else
            end
            signalDataSel{c} = signalData{condsToUse(c)}(repmat(selRowIx,[2,1]),:,curTrials); % repmat because the first half is real, second half is imag with same ordering
            noise1Sel{c}     =     noise1{condsToUse(c)}(repmat(selRowIx,[2,1]),:,curTrials);
            noise2Sel{c}     =     noise2{condsToUse(c)}(repmat(selRowIx,[2,1]),:,curTrials);
            
            % find non-empty frequency indices
            nonEmpty = find(cell2mat(cellfun(@(x) ~isempty(x),indF,'uni',false)));
            % find first among conditions to use
            nonEmpty = min(nonEmpty(ismember(nonEmpty,condsToUse)));
            % check if indices are unequal
            if any ( indF{nonEmpty} ~= indF{condsToUse(c)} )
                error('frequency indices are not matched across conditions');
            elseif any ( indB{nonEmpty} ~= indB{condsToUse(c)} )
                error('bin indices are not matched across conditions');
            else
            end

            % assign indices and labels of selected data features
            % grab bin indices
            indBSel{c} = indB{condsToUse(c)}(selRowIx);
            % grab frequency indices
            indFSel{c}  = indF{condsToUse(c)}(selRowIx);
            % grab bin labels
            binLabelsSel{c}  = binLabels{condsToUse(c)}(binsToUse+1); % add one, because bin level 0 = average
            % grap frequency labels
            freqLabelsSel{c}  = freqLabels{condsToUse(c)}(freqsToUse);
            % grap frequency labels
            trialsSel{c}  = curTrials;
        end
    end

    % now replace missing conditions with NaNs, dude!
    signalDataSel = replaceEmpty(signalDataSel);
    noise1Sel = replaceEmpty(noise1Sel);
    noise2Sel = replaceEmpty(noise2Sel);
    indFSel = replaceEmpty(indFSel);
    indBSel = replaceEmpty(indBSel);
    freqLabelsSel = replaceEmpty(freqLabelsSel);
    binLabelsSel = replaceEmpty(binLabelsSel);
    trialsSel = replaceEmpty(trialsSel);
end

function cellOut = replaceEmpty(cellIn)
    nonEmpty = cell2mat(cellfun(@(x) ~isempty(x),cellIn,'uni',false));
    repIdx = find(nonEmpty,1,'first');
    cellOut = cellIn;
    cellOut(~nonEmpty) = {nan(size(cellIn{repIdx}))};
end