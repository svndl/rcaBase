function [dataCube] = aggregateData(rcaData,freqIndices,keepConditions)
% [dataCube] = aggregateBins(rcaData,freqsToUse,[keepConditions])
%
% rcaData is any variable organized like the output of rcaRun ###
% freqIndices ###
%
% dataCube is a struct containing subject and trial averaged data that
% contains the following fields:
%
%   realBins: bin-by-harmonic-by-component array of REAL RC coefficients
%   imagBins: bin-by-harmonic-by-component array of IMAGINARY RC coefficients
%   ampBins: bin-by-harmonic-by-component array of RC amplitudes
%   phaseBins: bin-by-harmonic-by-component array of RC phases
%
%
% if keepConditions is true (default = false), condition number is the 4th
% dimension of each array.
%%
if nargin<3, keepConditions=false; end

%% convert to real/imaginary
[rcaDataReal,rcaDataImag] = getRealImag(rcaData);

if ~keepConditions
    % aggregate across subjects AND conditions
    rcaDataRealAll=cat(3,rcaDataReal{:});
    rcaDataImagAll=cat(3,rcaDataImag{:});
    
    % average over trials
    muRcaDataReal=nanmean(rcaDataRealAll,3);
    muRcaDataImag=nanmean(rcaDataImagAll,3);
    
    % make a bin by frequency by component cube
    nFreqs=numel(unique(freqIndices));
    nBins=size(muRcaDataReal,1)/nFreqs;
    realBins=zeros(nBins,nFreqs,nComp);
    imagBins=zeros(nBins,nFreqs,nComp);
    for c=1:nComp
        for f=1:nFreqs
            realBins(:,f,c)=muRcaDataReal(freqIndices==freqsToUse(f),c);
            imagBins(:,f,c)=muRcaDataImag(freqIndices==freqsToUse(f),c);
        end
    end

end

ampBins=sqrt(realBins.^2+imagBins.^2);
phaseBins=atan(imagBins./realBins);

dataCube.realBins=realBins;
dataCube.imagBins=imagBins;
dataCube.ampBins=ampBins;
dataCube.phaseBins=phaseBins;

% snrBins=2*ampBins./(noise1AmpBins+noise2AmpBins);
% ozSnrBins=2*ozAmpBins./(ozNoise1AmpBins+ozNoise2AmpBins);

