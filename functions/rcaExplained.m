function [ rcaRelExpl, rcaVarExpl, pcaVarExpl, pcaEigs ] = rcaExplained(rcStruct,nInclude,plotFig)
    if nargin < 3
        plotFig = false;
    else
    end
    if nargin < 2
        nInclude = 1;
    else
    end
    
    % compute PCA variance explained
    [~,pcaEigs] = eig(rcStruct.covData.Rxx + rcStruct.covData.Ryy);
    pcaEigs=sort(diag(pcaEigs),'descend');
    pcaVarExpl = sum(pcaEigs(1:nInclude))/sum(pcaEigs);
        
    % and reliability explained
    rcaEigs = sort(rcStruct.covData.sortedGeneralizedEigenValues,'descend');
    rcaRelExpl = sum(rcaEigs(1:nInclude))/sum(rcaEigs);

    % compute RCA variance explained
    Rtotal = 0.5*(rcStruct.covData.Rxx + rcStruct.covData.Ryy); % pooled covariance 128 x 128
    sigmas = diag (rcStruct.W'*Rtotal*rcStruct.W)  ./ diag( rcStruct.W'*rcStruct.W) ; % this should yield C x 1 matrix, where C is nComp 
    rcaVarExpl = sum(sigmas(1:nInclude))/sum(sigmas);
    if plotFig
        figure; plot(pcaEigs,'-o');
        xlim([0.5,50.5]); % only show first 50 components
    else 
    end   
end