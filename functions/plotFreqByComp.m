function [figNums] = plotFreqByComp(mainRcaData,mainNoiseData,rcaSettings,plotSettings,comparisonRcaData,comparisonNoiseData)
% [figNums] = plotFreqByComp(mainRcaData,mainNoiseData,rcaSettings,[plotSettings],[comparisonRcaData],[comparisonNoiseData])
%
% Create a numberFrequencies x numberRCs multipanel figure where the
% amplitude in microVolts is plotted against the bin levels. Optionally
% include an extra column for the comparison channel's data (e.g. for Oz).
% 
% Individual plots are modelled after PowerDiva with unfilled squares to
% indicate noise estimates for each bin.
%
% If plotSettings includes parameters for the conditions, then each
% condition will be plotted separately.


if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if (nargin<4 || isempty(plotSettings)), useSpecialSettings = false; else useSpecialSettings = true; end
if nargin<6, plotComparison = false; else plotComparison = true; end

poolOverBins = false;

separateConds = false;
allCondsInd = 1;
nCond = length(rcaSettings.condsToUse);
if useSpecialSettings
    if plotSettings.showConditions
        separateConds = true;
        allCondsInd = 1:nCond;
    else
        plotSettings.conditionColors = [0 0 0];
    end
    errorType = plotSettings.errorType;
else
    errorType = [];
    plotSettings.conditionColors = [0 0 0];
end

figNums = [];

nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);
nBinsHalf = floor(nBins/2);
binsToTick = [1 nBinsHalf nBins];

avgRcaData = aggregateData(mainRcaData,rcaSettings,separateConds,errorType);
avgNoise1Data = aggregateData(mainNoiseData.lowerSideBand,rcaSettings,separateConds,errorType);
avgNoise2Data = aggregateData(mainNoiseData.higherSideBand,rcaSettings,separateConds,errorType);

[~,noiseLevsMain] = computeSnr(avgRcaData,avgNoise1Data,avgNoise2Data,poolOverBins);

if plotComparison
    avgCompData = aggregateData(comparisonRcaData,rcaSettings,separateConds,errorType);
    avgCompNoise1Data = aggregateData(comparisonNoiseData.lowerSideBand,rcaSettings,separateConds,errorType);
    avgCompNoise2Data = aggregateData(comparisonNoiseData.higherSideBand,rcaSettings,separateConds,errorType);
    
    [~,noiseLevsCompare] = computeSnr(avgCompData,avgCompNoise1Data,avgCompNoise2Data,poolOverBins);
end

figure
for f=1:nFreqs
    for rc=1:rcaSettings.nComp
        
        if plotComparison
            subplot(nFreqs,rcaSettings.nComp+1,(f-1)*(rcaSettings.nComp+1)+rc);
        else
            subplot(nFreqs,rcaSettings.nComp,(f-1)*(rcaSettings.nComp)+rc);
        end
        
        for condNum = allCondsInd            
            plot(avgRcaData.ampBins(:,f,rc,condNum),'k-','LineWidth',1.5,'Color',plotSettings.conditionColors(condNum,:));
            hold on
            plot(noiseLevsMain(:,f,rc,condNum),'ks','Color',plotSettings.conditionColors(condNum,:));
            if ~isempty(errorType)                
                lb = avgRcaData.ampBins(:,f,rc,condNum) - avgRcaData.ampErrBins(:,f,rc,condNum,1);
                ub = avgRcaData.ampErrBins(:,f,rc,condNum,2) - avgRcaData.ampBins(:,f,rc,condNum);
                errorbar(1:nBins,avgRcaData.ampBins(:,f,rc,condNum),lb,ub,'Color',plotSettings.conditionColors(condNum,:));
            end
        end
        set(gca,'XTick',binsToTick);
        set(gca,'XTickLabel',round(rcaSettings.binLevels{condNum}(binsToTick)*100)./100); % ### assumes all conditions plotted have same binLevels
        xlim([0 nBins+1]);
        
        if f==1, title(['RC' num2str(rc)']); end
        
        if rc==1, ylabel(rcaSettings.freqLabels{f}); end
        
        if plotComparison
            subplot(nFreqs,rcaSettings.nComp+1,f*(rcaSettings.nComp+1));
            
            for condNum = allCondsInd
                plot(avgCompData.ampBins(:,f,1,condNum),'r-','LineWidth',1.5,'Color',plotSettings.conditionColors(condNum,:));
                hold on
                plot(noiseLevsCompare(:,f,1,condNum),'rs','Color',plotSettings.conditionColors(condNum,:));
                if ~isempty(errorType)
                    lb = avgCompData.ampBins(:,f,1,condNum) - avgCompData.ampErrBins(:,f,1,condNum,1);
                    ub = avgCompData.ampErrBins(:,f,1,condNum,2) - avgCompData.ampBins(:,f,1,condNum);
                    errorbar(1:nBins,avgCompData.ampBins(:,f,1,condNum),lb,ub,'Color',plotSettings.conditionColors(condNum,:));
                end
            end
            set(gca,'XTick',binsToTick);
            set(gca,'XTickLabel',round(rcaSettings.binLevels{condNum}(binsToTick)*100)./100); % ### assumes all conditions plotted have same binLevels
            xlim([0 nBins+1]);
            xlim([0 nBins+1]);
            if f==1
                if useSpecialSettings
                    title(plotSettings.comparisonName)
                else
                    title('Comparison Chan.')
                end
            end
        end
        
    end
end
figNums = [figNums,gcf];





