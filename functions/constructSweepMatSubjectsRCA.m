function [sweepMatSubjects] = constructSweepMatSubjectsRCA(rcaData,rcaSettings,noiseData1,noiseData2,compNum,condNum,freqNum)

nConditions = size(rcaData,1);
nSubjects = size(rcaData,2);
nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);
nCompFromInputData = size(rcaData{1,1},2);

sweepMatSubjects = nan(nBins,6,nSubjects);

signalData = getSubjectData(rcaData,rcaSettings,compNum,condNum,freqNum);

sweepMatSubjects(:,1:2,:) = signalData;

noiseData1 = getSubjectData(noiseData1,rcaSettings,compNum,condNum,freqNum);
sweepMatSubjects(:,3:4,:) = noiseData1;

noiseData2 = getSubjectData(noiseData2,rcaSettings,compNum,condNum,freqNum);
sweepMatSubjects(:,5:6,:) = noiseData2;
