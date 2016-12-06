function newH = grabCovFig(figureH)
    rcaH=get(figureH,'children');
    newH = figure;
    %ah=copyobj(rcaH(1),newH);
    figure(newH);
    axH=copyobj(rcaH(1),newH);
    %axH = findobj(rcaH(1),'type','axes');
    subplot(2,1,1,axH);
    axH=copyobj(rcaH(2),newH);
    subplot(2,1,2,axH);
    %set(ah,'position',[.1,.2,.8,.3]);
    %ah=copyobj(rcaH(2),newH);
    %set(ah,'position',[.1,.6,.8,.6]);
end

