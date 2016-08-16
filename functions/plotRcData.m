function [figNums,ampVals,noiseVals,binVals,threshInfo] = plotRcData(mainRcaData,mainNoiseData,rcaSettings,plotSettings,RCtoPlot,FreqToPlot,plotThreshold)
% [figNums] = plotRcData(mainRcaData,mainNoiseData,rcaSettings,plotSettings,RCtoPlot)
%
% Plot modelled after PowerDiva with squares to
% indicate noise estimates for each bin.
%
% If plotSettings includes parameters for the conditions, then each
% condition will be plotted separately.


if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if (nargin<4 || isempty(plotSettings)), useSpecialSettings = false; else useSpecialSettings = true; end
if nargin<5 || isempty(RCtoPlot), RCtoPlot = 1; end
if nargin<6 || isempty(FreqToPlot), FreqToPlot = 1; end
if nargin<7, plotThreshold = 0; end

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

nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);

if plotThreshold
    threshFitted = zeros(nFreqs,nCond);
    threshVal = nan(nFreqs,nCond);
    slopeVal = nan(nFreqs,nCond);
    fitBinRange = nan(nFreqs,nCond,2);
else
    threshVal = nan;
    slopeVal = nan;
    fitBinRange = nan;
end

avgRcaData = aggregateData(mainRcaData,rcaSettings,separateConds,errorType);
avgNoise1Data = aggregateData(mainNoiseData.lowerSideBand,rcaSettings,separateConds,errorType);
avgNoise2Data = aggregateData(mainNoiseData.higherSideBand,rcaSettings,separateConds,errorType);

[~,noiseLevsMain] = computeSnr(avgRcaData,avgNoise1Data,avgNoise2Data,poolOverBins);

figNums = [];
ampVals = nan(nCond,nBins,nFreqs);
noiseVals = nan(nCond,nBins,nFreqs);

rc = RCtoPlot;
f = FreqToPlot;

figure;
binVals = rcaSettings.binLevels{allCondsInd(1)}; % ### use first condition to be plotted.
if isLogSpaced(binVals)
    set(gca,'XScale','log');
end
hold on;
set(gca,'FontSize',24);
set(gca,'FontWeight','light');
set(gca,'Color','w');
set(gcf,'Color','w')
if ~isempty(plotSettings.titleToUse)
    title(plotSettings.titleToUse);
end
ylabel('Amplitude (\muVolts)')
xlabel(plotSettings.xlabel);

fThresh = fopen('ThresholdParameters.txt', 'a+');
for condNum = allCondsInd
    
    binVals = rcaSettings.binLevels{condNum};
    
    ampVals(condNum,:,f) = avgRcaData.ampBins(:,f,rc,condNum);
    noiseVals(condNum,:,f) = noiseLevsMain(:,f,rc,condNum);
    
    if isLogSpaced(binVals)
        semilogx(binVals,ampVals(condNum,:,f),'k-','LineWidth',3,'Color',plotSettings.conditionColors(condNum,:));
        semilogx(binVals,noiseVals(condNum,:,f),'ks','Color',plotSettings.conditionColors(condNum,:),'MarkerSize',12)%,...
        %'MarkerFaceColor',plotSettings.conditionColors(condNum,:));
    else
        plot(binVals,ampVals(condNum,:,f),'k-','LineWidth',3,'Color',plotSettings.conditionColors(condNum,:));
        plot(binVals,noiseVals(condNum,:,f),'ks','Color',plotSettings.conditionColors(condNum,:),'MarkerSize',12)%,...
        %'MarkerFaceColor',plotSettings.conditionColors(condNum,:));
    end
    if ~isempty(errorType)
        lb = ampVals(condNum,:,f)' - avgRcaData.ampErrBins(:,f,rc,condNum,1);
        ub = avgRcaData.ampErrBins(:,f,rc,condNum,2) - ampVals(condNum,:,f)';
        errorbar(binVals,ampVals(condNum,:,f),lb,ub,'Color',plotSettings.conditionColors(condNum,:),'LineWidth',1.5);
    end
    
    if plotThreshold
        clear sweepMatSubjects;
        sweepMatSubjects = constructSweepMatSubjectsRCA(mainRcaData,rcaSettings,...
            mainNoiseData.lowerSideBand,mainNoiseData.higherSideBand,rc,condNum,f);
        
        [tThr,tThrStdErr,tSlp,tSlpStdErr,tLSB,tRSB,~,tYFitPos,tXX, pValue] = getThreshScoringOutput(sweepMatSubjects, rcaSettings.binLevels{condNum});
        
        sPhases = sprintf(' %f ', avgRcaData.phaseBins(:, f, rc, condNum)');
        sAmps = sprintf(' %f ', ampVals(condNum,:,f)');        
        sBins = sprintf(' %f ', binVals);
        sAmpsLow = sprintf(' %f ', avgRcaData.ampErrBins(:, f, rc, condNum, 1)');
        sAmpsHigh = sprintf(' %f ', avgRcaData.ampErrBins(:, f, rc, condNum, 2)');
        pVals = sprintf(' %f ', pValue');
        fprintf(fThresh, '============================================ Experiment %s condition #%d  freq (%s) comp %s================================================\n', ...
            plotSettings.groupName, condNum,rcaSettings.freqLabels{f}, plotSettings.titleToUse);
        
        fprintf(fThresh, 'Bins:                  %s\n', sBins);
        fprintf(fThresh, 'Amplitude w err low:   %s\n', sAmpsLow);
        fprintf(fThresh, 'Amplitude w err high:  %s\n', sAmpsHigh);
        fprintf(fThresh, 'Phases:                %s\n', sPhases);
        fprintf(fThresh, 'Amplitude:             %s\n', sAmps);
        fprintf(fThresh, 'PValue:                %s\n', pVals);
        
        
        fprintf(fThresh, '\n');

        threshVal(f,condNum) = tThr;
        threshErr(f,condNum) = tThrStdErr;
        slopeVal(f,condNum) = tSlp;
        slopeErr(f,condNum) = tSlpStdErr;
        fitBinRange(f,condNum,:) = [tLSB,tRSB];
        if isnan(tThr)
            fprintf(fThresh, 'No threshold could be fitted for CondNum = %d (%s).\n',condNum,rcaSettings.freqLabels{f});
        else
            % save line info to plot after everything else so it's "on top"
            fprintf(fThresh, 'Thresh = %1.2f, Thresh Error = %1.2f, Slope = %1.2f, Slope Error = %1.2f, Range=[%d,%d], for CondNum = %d (%s).\n', ...
                threshVal(f,condNum), threshErr(f,condNum), ...
                slopeVal(f,condNum), slopeErr(f,condNum), ...
                fitBinRange(f,condNum,:), condNum, rcaSettings.freqLabels{f});
            threshFitted(f,condNum) = 1;
            saveXX{f,condNum} = tXX;
            saveY{f,condNum} = tYFitPos;
        end
        fprintf('\n');
    end
    
end

yoffset = [0.75 0.85];
if plotThreshold && any(threshFitted(f,:)>0)
    for condNum = 1:nCond
        binVals = rcaSettings.binLevels{condNum};
        if threshFitted(f,condNum)
            % add shaded error region on threshold values:
%             h = fill([threshVal(f,condNum)-threshErr(f,condNum) threshVal(f,condNum)+threshErr(f,condNum) threshVal(f,condNum)+threshErr(f,condNum) threshVal(f,condNum)-threshErr(f,condNum)],...
%                 [min(ylim) min(ylim) max(ylim) max(ylim)],'k');
%             set(h,'LineStyle','none','FaceAlpha',0.2,'FaceColor',plotSettings.conditionColors(condNum,:));
            
            yExtremes = [min(ylim) max(ylim)];
            plot(repmat(threshVal(f,condNum)-threshErr(f,condNum),[1 2]),...
                yExtremes,'k--','Color',plotSettings.conditionColors(condNum,:))
            plot(repmat(threshVal(f,condNum)+threshErr(f,condNum),[1 2]),...
                yExtremes,'k--','Color',plotSettings.conditionColors(condNum,:))
            
            % plot the linear fits & threshold values:
            if isLogSpaced(binVals)
                set(gca,'XScale','log');
                semilogx(saveXX{f,condNum},saveY{f,condNum},'k-','LineWidth',3);
                semilogx(threshVal(f,condNum),0,'kd','MarkerSize',18,...
                    'MarkerFaceColor',plotSettings.conditionColors(condNum,:),'LineWidth',3);
            else
                plot(saveXX{f,condNum},saveY{f,condNum},'k-','LineWidth',3);
                plot(threshVal(f,condNum),0,'kd','MarkerSize',18,...
                    'MarkerFaceColor',plotSettings.conditionColors(condNum,:),'LineWidth',3);
            end
            text(binVals(1),yoffset(condNum)*max(yExtremes),sprintf('thresh (slope): %2.3f+/-%2.3f (%2.3f+/-%2.3f)',threshVal(f,condNum),threshErr(f,condNum),slopeVal(f,condNum),slopeErr(f,condNum)),'Color',plotSettings.conditionColors(condNum,:),'FontSize',12);
            ylim(yExtremes);
            %text(.9*threshVal(f,condNum),.8*max(ylim),sprintf('%2.2f',threshVal(f,condNum)),'Color',plotSettings.conditionColors(condNum,:));
        end
    end
end

%set(gca,'XTick',plotSettings.xTick);
xlim([binVals(1)-(binVals(2)-binVals(1)) binVals(end)+(binVals(end)-binVals(end-1))])
%set(gca,'XTickLabel',round(rcaSettings.binLevels*10)./10);
if isfield(plotSettings,'ymax')
    ylim([0 plotSettings.ymax]);
end

% if plotThreshold && any(threshFitted(f,:)>0)
%     set(gca,'XTick',[min(xlim) binVals([1 floor(length(binVals)/2) end])])
% else
    set(gca,'XTick',binVals([1 floor(length(binVals)/2) end]))
% end

set(gca,'ticklength',1.75*get(gca,'ticklength'))

figNums = [figNums,gcf];

%text(0.8*max(xlim),0.8*max(ylim),rcaSettings.freqLabels{f},'FontSize',20);

if plotThreshold && any(threshFitted(f,:)>0)
    box on;
    threshInfo.threshFitted = threshFitted;
    threshInfo.xx = saveXX;
    threshInfo.YY = saveY;
    threshInfo.threshVals = threshVal;
    threshInfo.slopeVals = slopeVal;
    threshInfo.fitBinRange = fitBinRange;
end
fclose(fThresh);
