function [boo] = createSourceDataMat(dataPath)
%
% createSourceDataMat(dataPath)
%
% gather all relevant data from the PowerDiva exports & save it into a
% matlab variable on disk that contains data for all channels, all
% frequencies, all bins, all conditions
%
% datapath containing 'RLS_c*.txt' or 'DFT_c*.txt' files required

possDataTypes = {'RLS' 'DFT'};
boo = false;
for dT = 1:length(possDataTypes)
    
    dataType = possDataTypes{dT};
    filenames = dir([dataPath,'/',sprintf('%s*.txt',dataType)]);
    if ~isempty(filenames)
        fprintf('Reading in all %s data in %s ...  ',dataType,dataPath);                
        
        [signalData,indF,indB,noise1,noise2,freqLabels,binLevels,chanIncluded]=textExportToRca(dataPath,dataType);        
        
        sourceDataFileName = sprintf('%s/sourceData_%s.mat',dataPath,dataType);
        save(sourceDataFileName,'signalData','indF','indB','noise1','noise2','freqLabels','binLevels','chanIncluded');
        clear signalData indF indB noise1 noise2 freqLabels binLevels chanIncluded
        fprintf('Done.\n')
    else
        warning('Could not run, data not found: %s',[dataPath,'/',sprintf('%s*.txt',dataType)]);
    end

end
boo = true;