function [figNums] = plotCondByComp(mainRcaData,mainNoiseData,rcaSettings,plotSettings,comparisonRcaData,comparisonNoiseData,oneRC,showLegend)
% [figNums] = plotFreqByComp(mainRcaData,mainNoiseData,rcaSettings,[plotSettings],[comparisonRcaData],[comparisonNoiseData],[oneRC])
%
% Create a numberConditions x numberRCs multipanel figure where the
% amplitude in microVolts is plotted against the bin levels. Optionally
% include an extra column for the comparison channel's data (e.g. for Oz).
% 
% Individual plots are modelled after PowerDiva with unfilled squares to
% indicate noise estimates for each bin.
%


if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if (nargin<4 || isempty(plotSettings)), useSpecialSettings = false; else useSpecialSettings = true; end
if nargin<5 || isempty(comparisonRcaData), plotComparison = false; else plotComparison = true; end
if nargin<7 || isempty(oneRC), oneRC = 0; end

poolOverBins = false;
nCond = length(rcaSettings.condsToUse);
separateConds = true;
allCondsInd = 1:nCond;
errorType = plotSettings.errorType;

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

if ~oneRC
    rcsToPlot = 1:rcaSettings.nComp;
    lw = 1.5;
else
    rcsToPlot = oneRC;
    lw = 3;
end

figure;
set(gcf,'Color','w')
set(gca,'FontSize',18);
hold on;
for condNum = allCondsInd  
    for rc=rcsToPlot
        
        if ~oneRC
            if plotComparison
                subplot(nCond,rcaSettings.nComp+1,(condNum-1)*(rcaSettings.nComp+1)+rc);
            else
                subplot(nCond,rcaSettings.nComp,(condNum-1)*(rcaSettings.nComp)+rc);
            end
        else
            subplot(1,nCond,condNum);            
            set(gca,'FontSize',18);
            hold on;
            if condNum ==1
                ylabel('Amplitude (\muVolts)');
            end
        end
        
        for f=1:nFreqs    
            sc = 1-f/nFreqs;
            if ~mod(f,2)
                plot(avgRcaData.ampBins(:,f,rc,condNum),'k--','LineWidth',lw,'Color',sc.*plotSettings.conditionColors(condNum,:));
            else
                plot(avgRcaData.ampBins(:,f,rc,condNum),'k-','LineWidth',lw,'Color',sc.*plotSettings.conditionColors(condNum,:));
            end
            hold on
        end
        if ~oneRC
            if (rc==3), legend(rcaSettings.freqLabels,'Location','NorthWest'); end
        end
        for f=1:nFreqs
            sc = 1-f/nFreqs;
            %plot(noiseLevsMain(:,f,rc,condNum),'ks','Color',sc.*plotSettings.conditionColors(condNum,:));
            if ~isempty(errorType)                
                lb = avgRcaData.ampBins(:,f,rc,condNum) - avgRcaData.ampErrBins(:,f,rc,condNum,1);
                ub = avgRcaData.ampErrBins(:,f,rc,condNum,2) - avgRcaData.ampBins(:,f,rc,condNum);
                for b = 1:nBins
%                     if ~mod(f,2)
%                         errorbar(b,avgRcaData.ampBins(b,f,rc,condNum),lb(b),ub(b),...
%                             'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1);
%                     else
%                         errorbar(b,avgRcaData.ampBins(b,f,rc,condNum),lb(b),ub(b),...
%                             'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1.75);
%                     end                        
                end
            end
        end
        set(gca,'XTick',binsToTick);
        set(gca,'XTickLabel',round(rcaSettings.binLevels{condNum}(binsToTick)*100)./100); % ### assumes all conditions plotted have same binLevels
        xlim([0 nBins+1]);
        
        if ~oneRC
            if rc==1, ylabel(['Condition ',num2str(condNum),' Ampl.']); end
            if condNum==1, title(['RC' num2str(rc)']); end
        end
        
        if plotComparison && ~oneRC
            subplot(nCond,rcaSettings.nComp+1,condNum*(rcaSettings.nComp+1));
            
            for f=1:nFreqs
                sc = 1-f/nFreqs;
                if ~mod(f,2)
                    plot(avgCompData.ampBins(:,f,1,condNum),'r--','LineWidth',lw,'Color',sc.*plotSettings.conditionColors(condNum,:));
                else
                    plot(avgCompData.ampBins(:,f,1,condNum),'r-','LineWidth',lw,'Color',sc.*plotSettings.conditionColors(condNum,:));
                end
                hold on
                %plot(noiseLevsCompare(:,f,1,condNum),'rs','Color',sc.*plotSettings.conditionColors(condNum,:));
                if ~isempty(errorType)
                    lb = avgCompData.ampBins(:,f,1,condNum) - avgCompData.ampErrBins(:,f,1,condNum,1);
                    ub = avgCompData.ampErrBins(:,f,1,condNum,2) - avgCompData.ampBins(:,f,1,condNum);
                    for b = 1:nBins
%                         if ~mod(f,2)
%                             errorbar(b,avgCompData.ampBins(b,f,1,condNum),lb(b),ub(b),...
%                                 'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1);
%                         else
%                             errorbar(b,avgCompData.ampBins(b,f,1,condNum),lb(b),ub(b),...
%                                 'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1.75);                            
%                         end
                    end
                end
            end
            set(gca,'XTick',binsToTick);
            set(gca,'XTickLabel',round(rcaSettings.binLevels{condNum}(binsToTick)*100)./100); % ### assumes all conditions plotted have same binLevels
            xlim([0 nBins+1]);
            xlim([0 nBins+1]);
            if condNum==1
                if useSpecialSettings
                    title(plotSettings.comparisonName)
                else
                    title('Comparison Chan.')
                end
            end
        end
        
    end
end

if oneRC && showLegend
    for condNum= allCondsInd
        subplot(1,nCond,condNum);
        legend(rcaSettings.freqLabels,'Location','NorthWest');
    end
end

% fix y-limits of each subplot
ymax = [];
for condNum = allCondsInd  
    for rc=rcsToPlot
        
        if ~oneRC
            if plotComparison
                subplot(nCond,rcaSettings.nComp+1,(condNum-1)*(rcaSettings.nComp+1)+rc);
            else
                subplot(nCond,rcaSettings.nComp,(condNum-1)*(rcaSettings.nComp)+rc);
            end
        else
            subplot(1,nCond,condNum);
        end
        ymax = [ymax max(ylim)];
    end
end
globalMax = max(ymax);
for condNum = allCondsInd  
    for rc=rcsToPlot
        
        if ~oneRC
            if plotComparison
                subplot(nCond,rcaSettings.nComp+1,(condNum-1)*(rcaSettings.nComp+1)+rc);
            else
                subplot(nCond,rcaSettings.nComp,(condNum-1)*(rcaSettings.nComp)+rc);
            end
        else
            subplot(1,nCond,condNum);
        end
        ylim([0 globalMax]);
    end
end


figNums = [figNums,gcf];





