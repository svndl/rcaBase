function [dataCube] = aggregateOverSubjects(rcaData,freqsToUse)

% bindata is a struct containing subject and condition averaged data that
% contains the following fields:
%
%   realBins: bin-by-harmonic-by-component array of REAL RC coefficients
%   imagBins: bin-by-harmonic-by-component array of IMAGINARY RC coefficients
%   ampBins: bin-by-harmonic-by-component array of RC amplitudes
%   phaseBins: bin-by-harmonic-by-component array of RC phases
%
%   noise1RealBins: bin-by-harmonic-by-component array of noise-band (lo) real RC coefficients
%   noise1ImagBins: bin-by-harmonic-by-component array of noise-band (lo) imaginary RC coefficients
%   noise1AmpBins: bin-by-harmonic-by-component array of noise-band (lo) RC amplitudes
%   noise1PhaseBins: bin-by-harmonic-by-component array of noise-band (lo) RC phases
%
%   noise2RealBins: bin-by-harmonic-by-component array of noise-band (hi) real RC coefficients
%   noise2ImagBins: bin-by-harmonic-by-component array of noise-band (hi) imaginary RC coefficients
%   noise2AmpBins: bin-by-harmonic-by-component array of noise-band (hi) RC amplitudes
%   noise2PhaseBins: bin-by-harmonic-by-component array of noise-band (hi) RC phases
%

%% convert to real/imaginary
[rcaDataReal,rcaDataImag] = getRealImag(rcaData);
[noiseData1Real,noiseData1Imag] = getRealImag(noiseData.lowerSideBand);
[noiseData2Real,noiseData2Imag] = getRealImag(noiseData.higherSideBand);

%% aggregate across subjects AND conditions
rcaDataRealAll=cat(3,rcaDataReal{:});
rcaDataImagAll=cat(3,rcaDataImag{:});

noiseData1RealAll=cat(3,noiseData1Real{:});
noiseData1ImagAll=cat(3,noiseData1Imag{:});

noiseData2RealAll=cat(3,noiseData2Real{:});
noiseData2ImagAll=cat(3,noiseData2Imag{:});

ozDataRealAll=cat(3,ozDataReal{:});
ozDataImagAll=cat(3,ozDataImag{:});

ozNoiseData1RealAll=cat(3,ozNoiseData1Real{:});
ozNoiseData1ImagAll=cat(3,ozNoiseData1Imag{:});

ozNoiseData2RealAll=cat(3,ozNoiseData2Real{:});
ozNoiseData2ImagAll=cat(3,ozNoiseData2Imag{:});

%% trial average for sweep plots
muRcaDataReal=nanmean(rcaDataRealAll,3);
muRcaDataImag=nanmean(rcaDataImagAll,3);

muNoiseData1Real=nanmean(noiseData1RealAll,3);
muNoiseData1Imag=nanmean(noiseData1ImagAll,3);

muNoiseData2Real=nanmean(noiseData2RealAll,3);
muNoiseData2Imag=nanmean(noiseData2ImagAll,3);

muOzDataReal=nanmean(ozDataRealAll,3);
muOzDataImag=nanmean(ozDataImagAll,3);

muOzNoiseData1Real=nanmean(ozNoiseData1RealAll,3);
muOzNoiseData1Imag=nanmean(ozNoiseData1ImagAll,3);

muOzNoiseData2Real=nanmean(ozNoiseData2RealAll,3);
muOzNoiseData2Imag=nanmean(ozNoiseData2ImagAll,3);

%% make a frequency by bin by trial cube
nFreqs=numel(unique(freqIndices));
nBins=numel(unique(binIndices));
realBins=zeros(nBins,nFreqs,nComp);
imagBins=zeros(nBins,nFreqs,nComp);
noise1RealBins=zeros(nBins,nFreqs,nComp);
noise1ImagBins=zeros(nBins,nFreqs,nComp);
noise2RealBins=zeros(nBins,nFreqs,nComp);
noise2ImagBins=zeros(nBins,nFreqs,nComp);
ozRealBins=zeros(nBins,nFreqs,nComp);
ozImagBins=zeros(nBins,nFreqs,nComp);
ozNoise1RealBins=zeros(nBins,nFreqs,nComp);
ozNoise1ImagBins=zeros(nBins,nFreqs,nComp);
ozNoise2RealBins=zeros(nBins,nFreqs,nComp);
ozNoise2ImagBins=zeros(nBins,nFreqs,nComp);
for c=1:nComp
    for f=1:nFreqs
        realBins(:,f,c)=muRcaDataReal(freqIndices==freqsToUse(f),c);
        imagBins(:,f,c)=muRcaDataImag(freqIndices==freqsToUse(f),c);
        
        noise1RealBins(:,f,c)=muNoiseData1Real(freqIndices==freqsToUse(f),c);
        noise1ImagBins(:,f,c)=muNoiseData1Imag(freqIndices==freqsToUse(f),c);
        
        noise2RealBins(:,f,c)=muNoiseData2Real(freqIndices==freqsToUse(f),c);
        noise2ImagBins(:,f,c)=muNoiseData2Imag(freqIndices==freqsToUse(f),c);
        
        ozRealBins(:,f,c)=muOzDataReal(freqIndices==freqsToUse(f),1); % just repeat for Oz
        ozImagBins(:,f,c)=muOzDataImag(freqIndices==freqsToUse(f),1);
        
        ozNoise1RealBins(:,f,c)=muOzNoiseData1Real(freqIndices==freqsToUse(f),1);
        ozNoise1ImagBins(:,f,c)=muOzNoiseData1Imag(freqIndices==freqsToUse(f),1);
        
        ozNoise2RealBins(:,f,c)=muOzNoiseData2Real(freqIndices==freqsToUse(f),1);
        ozNoise2ImagBins(:,f,c)=muOzNoiseData2Imag(freqIndices==freqsToUse(f),1);
        
    end
end

%% extra info

ampBins=sqrt(realBins.^2+imagBins.^2);
phaseBins=atan(imagBins./realBins);

noise1AmpBins=sqrt(noise1RealBins.^2+noise1ImagBins.^2);
noise1PhaseBins=atan(noise1ImagBins./noise1RealBins);

noise2AmpBins=sqrt(noise2RealBins.^2+noise2ImagBins.^2);
noise2PhaseBins=atan(noise2ImagBins./noise2RealBins);

ozAmpBins=sqrt(ozRealBins.^2+ozImagBins.^2);
ozPhaseBins=atan(ozImagBins./ozRealBins);

ozNoise1AmpBins=sqrt(ozNoise1RealBins.^2+ozNoise1ImagBins.^2);
ozNoise1PhaseBins=atan(ozNoise1ImagBins./ozNoise1RealBins);

ozNoise2AmpBins=sqrt(ozNoise2RealBins.^2+ozNoise2ImagBins.^2);
ozNoise2PhaseBins=atan(ozNoise2ImagBins./ozNoise2RealBins);

snrBins=2*ampBins./(noise1AmpBins+noise2AmpBins);
ozSnrBins=2*ozAmpBins./(ozNoise1AmpBins+ozNoise2AmpBins);

%%

binData.realBins=realBins;
binData.imagBins=imagBins;
binData.ampBins=ampBins;
binData.phaseBins=phaseBins;

binData.noise1RealBins=noise1RealBins;
binData.noise1ImagBins=noise1ImagBins;
binData.noise1AmpBins=noise1AmpBins;
binData.noise1PhaseBins=noise1PhaseBins;

binData.noise2RealBins=noise2RealBins;
binData.noise2ImagBins=noise2ImagBins;
binData.noise2AmpBins=noise2AmpBins;
binData.noise2PhaseBins=noise2PhaseBins;

binData.ozRealBins=ozRealBins;
binData.ozImagBins=ozImagBins;
binData.ozAmpBins=ozAmpBins;
binData.ozPhaseBins=ozPhaseBins;

binData.ozNoise1RealBins=ozNoise1RealBins;
binData.ozNoise1ImagBins=ozNoise1ImagBins;
binData.ozNoise1AmpBins=ozNoise1AmpBins;
binData.ozNoise1PhaseBins=ozNoise1PhaseBins;

binData.ozNoise2RealBins=ozNoise2RealBins;
binData.ozNoise2ImagBins=ozNoise2ImagBins;
binData.ozNoise2AmpBins=ozNoise2AmpBins;
binData.ozNoise2PhaseBins=ozNoise2PhaseBins;

binData.snrBins=snrBins;
binData.ozSnrBins=ozSnrBins;
