function [figNums] = plotSnr(mainRcaData,mainNoiseData,rcaSettings,plotSettings,comparisonRcaData,comparisonNoiseData)

% ### add functionality for condition separation

if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if nargin<4, useSpecialSettings = false; end
if (nargin>=4 && ~isempty(plotSettings)), useSpecialSettings = true; else useSpecialSettings = false; end
if nargin<6, plotComparison = false; else plotComparison = true; end

poolOverBins = true;

figNums = [];

nFreqs = length(rcaSettings.freqsToUse);

avgRcaData = aggregateData(mainRcaData,rcaSettings);
avgNoise1Data = aggregateData(mainNoiseData.lowerSideBand,rcaSettings);
avgNoise2Data = aggregateData(mainNoiseData.higherSideBand,rcaSettings);

snrMain = computeSnr(avgRcaData,avgNoise1Data,avgNoise2Data,poolOverBins);

if plotComparison
    avgCompData = aggregateData(comparisonRcaData,rcaSettings);
    avgCompNoise1Data = aggregateData(comparisonNoiseData.lowerSideBand,rcaSettings);
    avgCompNoise2Data = aggregateData(comparisonNoiseData.higherSideBand,rcaSettings);
    
    snrComparison = computeSnr(avgCompData,avgCompNoise1Data,avgCompNoise2Data,poolOverBins);
end

xticks = 1:length(rcaSettings.binLevels{1}); 
if plotComparison
    % plot SNR of RC1 and comparison channel, using noise estimates pooled over all bins
    figure;
    set(gca,'Color','w');
    for rc = 1:rcaSettings.nComp
        for f=1:nFreqs
            subplot(rcaSettings.nComp,nFreqs,f+nFreqs*(rc-1)); hold on 
            color = (rc/rcaSettings.nComp).*[.7 .7 .7];
            plot(snrMain(:,f,rc),'-ok','MarkerFaceColor',color,'MarkerEdgeColor','none','Color',color);
            dataLabels = sprintf('RC%d',rc);
            plot(snrComparison(:,f,1),'-or','MarkerFaceColor','r','MarkerEdgeColor','none');
            if useSpecialSettings
                dataLabels = {dataLabels,plotSettings.comparisonName};
            else
                dataLabels = {dataLabels,'Comparison'};
            end
            title(rcaSettings.freqLabels{f});
            if f==1
                ylabel('SNR');
            end
            set(gca,'XTickLabel',xticks([1 floor(length(xticks)/2) end]));
            if f==1, hlg=legend(dataLabels,'Location','NorthWest'); set(hlg,'box','off'); end
            axis square;
        end
    end
    figNums = [figNums,gcf];
end

% plot the SNRs of each RC computed, using noise estimates pooled over all bins
% organize subplots by frequency
figure;
set(gca,'Color','w');
for f=1:nFreqs
    for rc=1:rcaSettings.nComp
        subplot(1,nFreqs,f);
        color = (rc/rcaSettings.nComp).*[.7 .7 .7];
        plot(snrMain(:,f,rc),'-o','Color',color,'MarkerFaceColor',color);
        hold on;
        dataLabels{rc} = sprintf('RC%d',rc);
        title(rcaSettings.freqLabels{f});
        if f==1
            ylabel('SNR');
        end
        set(gca,'XTickLabel',xticks([1 floor(length(xticks)/2) end])); 
        axis square;
    end
    if f==1, hlg=legend(dataLabels,'Location','NorthWest'); end
    set(hlg,'box','off');
end
figNums = [figNums,gcf];


% plot the SNRs of each RC computed, using noise estimates pooled over all bins
% organize subplots by RC
figure;
set(gca,'Color','w');
for rc=1:rcaSettings.nComp
    for f=1:nFreqs
        subplot(1,rcaSettings.nComp,rc);
        color = (f/nFreqs).*[.7 .7 .7];
        plot(snrMain(:,f,rc),'-o','Color',color,'MarkerFaceColor',color);
        hold on;
        dataLabels{f} = rcaSettings.freqLabels{f}; 
        title(sprintf('RC%d',rc));
        if rc==1
            ylabel('SNR');
        end
        set(gca,'XTickLabel',xticks([1 floor(length(xticks)/2) end])); 
        axis square;
    end
    if rc==1, hlg=legend(dataLabels,'Location','NorthWest'); end
    set(hlg,'box','off');
end
figNums = [figNums,gcf];





