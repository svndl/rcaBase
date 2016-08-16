function [signalDataSel,indFSel,indBSel,noise1Sel,noise2Sel,freqLabelsSel,binLevelsSel] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse)
%
% [signalData,indF,indB,noise1,noise2] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,condsToUse)
%
% Return only the desired subset of data contained in sourceDataFileName
% which contains only binsToUse, freqsToUse, condsToUse

load(sourceDataFileName); % should contain these variables: signalData,indF,indB,noise1,noise2,freqLabels,binLevels,chanIncluded

nR = length(condsToUse);
nC = size(signalData,2);
signalDataSel = cell(nR,nC);
noise1Sel = cell(nR,nC);
noise2Sel = cell(nR,nC);
binLevelsSel = cell(nR,1);
selRowIx = ismember(indB,binsToUse) & ismember(indF,freqsToUse);
for row = 1:nR
    for col = 1:nC
        signalDataSel{row,col} = signalData{condsToUse(row),col}(repmat(selRowIx,[2 1]),:,:); % repmat because the first half is real, second half is imag with same ordering
        noise1Sel{row,col}     =     noise1{condsToUse(row),col}(repmat(selRowIx,[2 1]),:,:);
        noise2Sel{row,col}     =     noise2{condsToUse(row),col}(repmat(selRowIx,[2 1]),:,:);
        
        binLevelsSel{row} = binLevels{condsToUse(row)}(binsToUse);
    end
end

indBSel = indB(selRowIx);
indFSel = indF(selRowIx);

for k = 1:length(freqsToUse)
    freqLabelsSel{k,1} = freqLabels{freqsToUse(k)};
end