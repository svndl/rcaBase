function [snr,noiseLevs] = computeSnr(mainData,noise1Data,noise2Data,poolOverBins)

if nargin<4, poolOverBins = true; end

sigLevs = mainData.ampBins;

if poolOverBins
    noiseLevs = (noise1Data.ampBins+noise2Data.ampBins)./2;
    pooledNoiseOverBins = mean(noiseLevs);
    noiseLevs = repmat(pooledNoiseOverBins,[size(noiseLevs,1) 1]);
else
    noiseLevs = (noise1Data.ampBins+noise2Data.ampBins)./2;
end

snr = sigLevs./noiseLevs;