function [dataOut,W,A,C_s,C_n,dGen,plotSettings] = ssdRun_rcaFormat(data,noise1,noise2,nReg,nComp,condRange,subjRange,show,plotStyle)
% [DATAOUT,W,A,C_s,C_n,dGen,plotSettings]=RCARUN(DATA,NOISE1,NOISE2,DATARCAFORMAT,[NREG],[NCOMP],[CONDRANGE],[SUBJRANGE],[SHOW],[PLOTSTYLE],)
% perform SSD dimensionality reduction: learn SNR-maximizing filter
%   and project data into the corresponding space
% 
% data: data can either be a 3-D array (samples x channels x trials), or a cell
%   array (conditions x subjects) of 3-D data arrays (as above). Here,
%   samples corresponds to (sweep) bins * frequencies and each entry is a
%   complex number corresponding to the (DFT) coefficient of a frequency.
% noise1: lowerSideBandNoise, same dimensions as data.
% noise2: higherSideBandNoise, same dimensions as data.
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
% function takes data as for rca (real and imaginary part concatenated in one matrix).
% ssd is calculated on complex data, change of representation is done in this function
%
% file is adapted from rcaRun.m by Jacek P. Dmochowski, 2014
% SB, 2016

%% TODO:
% (1) check subplot numbering
% (2) check preComputeRcaCovariancesLoop for mean centering consistency

if nargin<9 || isempty(plotStyle), plotStyle = 'matchMaxSignsToRc1'; end;
if nargin<8 || isempty(show), show=1; end;
if nargin<7 || isempty(subjRange) 
    if iscell(data)
        subjRange=1:size(data,2);
    else
        subjRange=1;
    end
end
if nargin<6 || isempty(condRange) 
    if iscell(data)
        condRange=1:size(data,1);
    else
        condRange=1;
    end
end
if nargin<5 || isempty(nComp), nComp=3; end
if nargin<4 || isempty(nReg), nReg=[]; end
if nargin<3 || isempty(data) || isempty(noise1) || isempty(noise2), error('At least three arguments is required (signal and 2 * noise data)'); end

if nComp>nReg, error('number of components cannot exceed the regularization parameter'); end; 


% sensor data, noise1 and noise2 need to be converted to complex
% representation for training the ssd-filters
[dataReal,dataImag] =getRealImag(data) ;
dataCmplx            = getCmplx(dataReal, dataImag);

[noise1Real,noise1Imag] =getRealImag(noise1) ;
noise1Cmplx            = getCmplx(noise1Real, noise1Imag);

[noise2Real,noise2Imag] =getRealImag(noise2) ;
noise2Cmplx            = getCmplx(noise2Real, noise2Imag);



% compute un-normalized covariances and keep track of number non-NaN data points
% in each case
if iscell(data)
    [sum_S,sum_N,nPointsIn_S,nPointsIn_N]=preComputeSsdCovariancesLoop(dataCmplx,noise1Cmplx,noise2Cmplx,condRange,subjRange);
else
    [sum_S,sum_N,nPointsIn_S,nPointsIn_N]=preComputeSsdCovariancesLoop(dataCmplx,noise1Cmplx,noise2Cmplx);
end
    
% accumulate covariances across selected conditions/subjects
fprintf('Accumulating covariance across selected conditions/subjects... \n');
switch ndims(sum_S)
    case 2  % channels x channels
        C_s=sum_S./nPointsIn_S; C_n=sum_N./nPointsIn_N; 
    case 3 % cond x channels x channels OR subjects x channels x channels
        if length(subjRange)==1
            C_s=squeeze(sum(sum_S,1))./squeeze(sum(nPointsIn_S,1));
            C_n=squeeze(sum(sum_N,1))./squeeze(sum(nPointsIn_N,1));
        elseif length(condRange)==1
            C_s=squeeze(sum(sum_S,1))./squeeze(sum(nPointsIn_S,1));
            C_n=squeeze(sum(sum_N,1))./squeeze(sum(nPointsIn_N,1));
        else
            error('something went wrong')
        end
    case 4  % cond x subj x channels x channels
        if length(subjRange)>1 && length(condRange)>1
            C_s=squeeze(sum(squeeze(sum(sum_S,2)),1))./squeeze(sum(squeeze(sum(nPointsIn_S,2)),1));
            C_n=squeeze(sum(squeeze(sum(sum_N,2)),1))./squeeze(sum(squeeze(sum(nPointsIn_N,2)),1));
        elseif length(condRange)==1
            C_s=squeeze(sum(sum_S,2)) ./ squeeze(sum(nPointsIn_S,2)) ;
            C_n=squeeze(sum(sum_N,2)) ./ squeeze(sum(nPointsIn_N,2)) ;
        elseif length(subjRange)==1
            C_s=squeeze(sum(squeeze(sum_S),1))./squeeze(sum(squeeze(nPointsIn_S),1));
            C_n=squeeze(sum(squeeze(sum_N),1))./squeeze(sum(squeeze(nPointsIn_N),1));
        else error('something went wrong');    
            
        end   
end

% train SSD spatial filters
[W,A,~,dGen,Kout] = ssdTrain(C_s,C_n,nComp,nReg);
%data = dataRcaFormat;

% project data into SSD space
if iscell(data)
    [nCond,nSubjects]=size(data);
    dataOut=cell(nCond,nSubjects);
    for subj=1:nSubjects
        for cond=1:nCond
            fprintf('Projecting into SSD space for subject %d/%d and condition %d/%d... \n',subj,nSubjects,cond,nCond);
            thisVolume=data{cond,subj};
            % projection as such is the same method, only the filters W have been calculated differently
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
            plotOnEgi(s(c).*A(:,c),colorbarLimits);
            title(['SSD C' num2str(c)]);
            axis off;
        end
    catch
        fprintf('call to plotOnEgi() failed. Plotting electrode values in a 1D style instead.\n');
        for c=1:nComp, subplot(3,nComp,c); plot(A(:,c),'*k-'); end
        title(['SSD C' num2str(c)]);
    end
    
    try
        for c=1:nComp
            subplot(3,nComp,c+nComp);
            shadedErrorBar([],s(c).*muData(:,c),semData(:,c),'k');
            title(['SSD C' num2str(c) ' time course']);
            axis tight;
        end
    catch
        fprintf('unable to plot ssd c means and sems. \n');
    end
    
    try
        [~,eigs]=eig(C_s);
        eigs=sort(diag(eigs),'ascend');
        subplot(325); hold on
        plot(eigs,'*k:'); 
        nEigs=length(eigs);
        plot(nEigs-Kout,eigs(nEigs-Kout),'*g');
        title('Signal Eigenvalue spectrum');
    catch
        fprintf('unable to plot within-trial covariance spectrum. \n');
    end
    
    subplot(326); plot(dGen,'*k:'); title('Generalized Eigenvalue spectrum');
    
end

