function [signalDataSel,indFSel,indBSel,noise1Sel,noise2Sel,freqLabelsSel,binLevelsSel,trialsSel] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse)
%
% [signalData,indF,indB,noise1,noise2,freqLabelsSel,binLevelsSel,trialsSel] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse,trialsToUse)
%
% Return only the desired subset of data contained in sourceDataFileName
% which contains only binsToUse, freqsToUse, condsToUse, trialsToUse
    load(sourceDataFileName); % should contain these variables: signalData,indF,indB,noise1,noise2,freqLabels,binLevels,chanIncluded

    if nargin < 5; trialsToUse = []; else end
    if nargin < 4 || isempty(condsToUse); condsToUse = 1:size(signalData,2); else end;
    if nargin < 3 || isempty(condsToUse); freqsToUse = unique(indF); else end;
    if nargin < 2 || isempty(condsToUse); 
        binsToUse = unique(indB); 
        binsToUse = binsToUse(binsToUse>0); % use all bins, but not the average bin
    else
    end

    nC = size(signalData,2);
    nR = length(condsToUse);
    signalDataSel = cell(nR,nC);
    noise1Sel = cell(nR,nC);
    noise2Sel = cell(nR,nC);
    binLevelsSel = cell(nR,1);
    selRowIx = ismember(indB,binsToUse) & ismember(indF,freqsToUse);
    for row = 1:nR
        for col = 1:nC
            if isempty(trialsToUse)
                trialsSel = 1:size(signalData{condsToUse(row),col},3); % use all trials
            else
                trialsSel = trialsToUse;
            end
            missingIdx = find(~ismember(trialsSel,1:size(signalData{condsToUse(row),col},3)));
            if ~isempty(missingIdx)
                error('Input trial indices "%s" is not among set of trials',num2str(trialsSel(missingIdx),'%d,'));
            else
            end
            signalDataSel{row,col} = signalData{condsToUse(row),col}(repmat(selRowIx,[2 1]),:,trialsSel); % repmat because the first half is real, second half is imag with same ordering
            noise1Sel{row,col}     =     noise1{condsToUse(row),col}(repmat(selRowIx,[2 1]),:,trialsSel);
            noise2Sel{row,col}     =     noise2{condsToUse(row),col}(repmat(selRowIx,[2 1]),:,trialsSel);
        
            binLevelsSel{row} = binLevels{condsToUse(row)}(binsToUse+1); % add one, because bin level 0 = average
        end
    end

    indBSel = indB(selRowIx);
    indFSel = indF(selRowIx);

    for k = 1:length(freqsToUse)
        freqLabelsSel{k,1} = freqLabels{freqsToUse(k)};
    end
end
