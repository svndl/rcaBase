function [dataOut,W,A,Rxx,Ryy,Rxy,dGen,plotSettings] = rcaRun(data,nReg,nComp,condRange,subjRange,show,plotStyle)
% [DATAOUT,W,A,RXX,RYY,RXY]=RCARUN(DATA,[NREG],[NCOMP],[CONDRANGE],[SUBJRANGE],[SHOW],[PLOTSTYLE])
% perform RCA dimensionality reduction: learn reliability-maximizing filter
%   and project data into the corresponding space
% 
% data: a 3-D array (samples x channels x trials) % data can either be a 3-D array (samples x channels x trials), or a cell
%   array (conditions x subjects) of 3-D data arrays (as above)
% nReg: regularization parameter controlling the number of bases to
%   diagonalize pooled autocovariance (must be positive integer less than the
%   number of channels).  Defaults to the number of bases diagonalizing 60%
%   of pooled covariance.  For EEG, do not set this to be the number of
%   electrodes.  Typical values range from 5-15.  
% nComp: number of components to return (defaults to 3)
% condRange: if data is a cell array, a subset of conditions on which to
%   learn RCA spatial filters (defaults to all conditions)
% subjRange: if data is a cell array, a subset of subjects on which to
%   learn RCA spatial filters (defaults to all subjects)
% show: 1 (default) to show learned components, 0 to not show
% plotStyle: 'matchMaxSignsToRc1' (default) or 'orig'
%   'matchMaxSignsToRc1' : the sign of the maximum magnitude value
%       in each RC is matched to the sign of RC1's maximum magnitude value
%       so that polarity colors are visually easier to compare. Each RC has
%       identical symmetric colorbars.
%   'orig' : each RC's topomap has its own autoscaled colorbar, and signs
%       are preserved from the output of rcaTrain
%
% dataOut: if data is a cell, this is a corresponding cell array
% (conditions x subjects) of dimensionality-reduced data volumes (samples x
% components x trials); if data is a 3D array, dataOut is the corresponding
% array (again, samples x components x trials)
% W: matrix of learned optimal projections (channels x nComp)
% A: matrix of resulting "forward models" -- the projection of the reliable
%   activity extracted by W onto the original electrode space (see Parra et
%   al., 2004, Neuroimage)
%
% (c) Jacek P. Dmochowski, 2014

%% TODO:
% (1) check subplot numbering
% (2) check preComputeRcaCovariancesLoop for mean centering consistency

if nargin<7 || isempty(plotStyle), plotStyle = 'matchMaxSignsToRc1'; end;
if nargin<6 || isempty(show), show=1; end;
if nargin<5 || isempty(subjRange) 
    if iscell(data)
        subjRange=1:size(data,2);
    else
        subjRange=1;
    end
end
if nargin<4 || isempty(condRange) 
    if iscell(data)
        condRange=1:size(data,1);
    else
        condRange=1;
    end
end
if nargin<3 || isempty(nComp), nComp=3; end
if nargin<2 || isempty(nReg), nReg=[]; end
if nargin<1 || isempty(data), error('At least one argument is required'); end


if nComp>nReg, error('number of components cannot exceed the regularization parameter'); end; 


% compute un-normalized covariances and keep track of number non-NaN data points
% in each case
if iscell(data)
    %[sumXX,sumYY,sumXY,nPointsInXX,nPointsInYY,nPointsInXY]=preComputeRcaCovariances(data,condRange,subjRange);
    [sumXX,sumYY,sumXY,nPointsInXX,nPointsInYY,nPointsInXY]=preComputeRcaCovariancesLoop(data,condRange,subjRange);
else
    %[sumXX,sumYY,sumXY,nPointsInXX,nPointsInYY,nPointsInXY]=preComputeRcaCovariances(data);
    [sumXX,sumYY,sumXY,nPointsInXX,nPointsInYY,nPointsInXY]=preComputeRcaCovariancesLoop(data);
end
    


% accumulate covariances across selected conditions/subjects
fprintf('Accumulating covariance across selected conditions/subjects... \n');
switch ndims(sumXX)
    case 2  % channels x channels
        Rxx=sumXX./nPointsInXX; Ryy=sumYY./nPointsInYY; Rxy=sumXY./nPointsInXY; 
    case 3 % cond x channels x channels OR subjects x channels x channels
        if length(subjRange)==1
            Rxx=squeeze(sum(sumXX,1))./squeeze(sum(nPointsInXX,1));
            Ryy=squeeze(sum(sumYY,1))./squeeze(sum(nPointsInYY,1));
            Rxy=squeeze(sum(sumXY,1))./squeeze(sum(nPointsInXY,1));
        elseif length(condRange)==1
            Rxx=squeeze(sum(sumXX,1))./squeeze(sum(nPointsInXX,1));
            Ryy=squeeze(sum(sumYY,1))./squeeze(sum(nPointsInYY,1));
            Rxy=squeeze(sum(sumXY,1))./squeeze(sum(nPointsInXY,1));
        else
            error('something went wrong')
        end
    case 4  % cond x subj x channels x channels
        if length(subjRange)>1 && length(condRange)>1
            Rxx=squeeze(sum(squeeze(sum(sumXX,2)),1))./squeeze(sum(squeeze(sum(nPointsInXX,2)),1));
            Ryy=squeeze(sum(squeeze(sum(sumYY,2)),1))./squeeze(sum(squeeze(sum(nPointsInYY,2)),1));
            Rxy=squeeze(sum(squeeze(sum(sumXY,2)),1))./squeeze(sum(squeeze(sum(nPointsInXY,2)),1));
        elseif length(condRange)==1
            Rxx=squeeze(sum(sumXX,2)) ./ squeeze(sum(nPointsInXX,2)) ;
            Ryy=squeeze(sum(sumYY,2)) ./ squeeze(sum(nPointsInYY,2)) ;
            Rxy=squeeze(sum(sumXY,2)) ./ squeeze(sum(nPointsInXY,2)) ;
        elseif length(subjRange)==1
            Rxx=squeeze(sum(squeeze(sumXX),1))./squeeze(sum(squeeze(nPointsInXX),1));
            Ryy=squeeze(sum(squeeze(sumYY),1))./squeeze(sum(squeeze(nPointsInYY),1));
            Rxy=squeeze(sum(squeeze(sumXY),1))./squeeze(sum(squeeze(nPointsInXY),1));
        else error('something went wrong');    
            
        end
        
end

% train RCA spatial filters
[W,A,~,dGen,Kout] = rcaTrain(Rxx,Ryy,Rxy,nReg,nComp);

% project data into RCA space
if iscell(data)
    [nCond,nSubjects]=size(data);
    dataOut=cell(nCond,nSubjects);
    for subj=1:nSubjects
        for cond=1:nCond
            fprintf('Projecting into RCA space for subject %d/%d and condition %d/%d... \n',subj,nSubjects,cond,nCond);
            thisVolume=data{cond,subj};
            dataOut{cond,subj}=rcaProject(thisVolume,W);
        end
    end
    try
        catData=cat(3,dataOut{:});
        muData=nanmean(catData,3);
        semData=nanstd(catData,[],3)/sqrt(size(catData,3));
    catch
        fprintf('could not compute means and sems \n');
    end
else
    dataOut = rcaProject(data,W);
    muData=nanmean(dataOut,3);
    semData=nanstd(dataOut,[],3)/sqrt(size(dataOut,3));
end

if show
    h=figure;
    
    if strcmp(plotStyle,'matchMaxSignsToRc1')
        symmetricColorbars = true;
        alignPolarityToRc1 = true;
    else % use original plotting conventions
        symmetricColorbars = false;
        alignPolarityToRc1 = false;        
    end
    
    if symmetricColorbars
        % for a consistent colorbar across RCs:
        colorbarLimits = [min(A(:)),max(A(:))];
        newExtreme = max(abs(colorbarLimits));
        colorbarLimits = [-newExtreme,newExtreme];
    else
        colorbarLimits = [];
    end
    if alignPolarityToRc1
        extremeVals = [min(A); max(A)];
        for rc = 1:nComp
            [~,f(rc)]=max(abs(extremeVals(:,rc)));
        end
        s = ones(1,nComp);
        if f(1)==1 % largest extreme value of RC1 is negative
            s(1) = -1; % flip the sign of RC1 so its largest extreme value is positive (yellow in parula)
            f(1) = 2; 
        end
        for rc = 2:nComp
            if f(rc)~=f(1) % if the poles containing the maximum corr coef are different
                s(rc) = -1; % we will flip the sign of that RC's time course & thus its corr coef values (in A) too
            end
        end
    else
        s = ones(1,nComp);
    end
    
    plotSettings.signFlips = s;
    plotSettings.colorbarLimits = colorbarLimits;
    
    try
        for c=1:nComp
            subplot(3,nComp,c);
            if ~isempty(which('mrC.plotOnEgi')) % check for mrC version of plotOnEgi
                mrC.plotOnEgi(s(c).*A(:,c),colorbarLimits);
            else
                plotOnEgi(s(c).*A(:,c),colorbarLimits);
            end
            title(['RC' num2str(c)]);
            axis off;
        end
    catch
        fprintf('call to plotOnEgi() failed. Plotting electrode values in a 1D style instead.\n');
        for c=1:nComp, subplot(3,nComp,c); plot(A(:,c),'*k-'); end
        title(['RC' num2str(c)]);
    end
    
    try
        for c=1:nComp
            subplot(3,nComp,c+nComp);
            shadedErrorBar([],s(c).*muData(:,c),semData(:,c),'k');
            title(['RC' num2str(c) ' time course']);
            axis tight;
        end
    catch
        fprintf('unable to plot rc means and sems. \n');
    end
    
    try
        [~,eigs]=eig(Rxx+Ryy);
        eigs=sort(diag(eigs),'ascend');
        subplot(325); hold on
        plot(eigs,'*k:'); 
        nEigs=length(eigs);
        plot(nEigs-Kout,eigs(nEigs-Kout),'*g');
        title('Within-trial covariance spectrum');
    catch
        fprintf('unable to plot within-trial covariance spectrum. \n');
    end
    
    subplot(326); plot(dGen,'*k:'); title('Across-trial covariance spectrum');
    
end

