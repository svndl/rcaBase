%% Example script to call rcaSweep on a dataset
% 
% Note! To run the full functionality of this code, you must also have the 
% free, public repository of code called sweepAnalysis/ which includes helper
% functions that are used for estimating error ellipses, confidence intervals,
% threshold fitting, jackknifing errors, and more. Please be sure to download 
% the latest version from:
% 
% https://github.com/hgerhard/sweepAnalysis
%
% Be sure to add the directory to your Matlab path, so that it can find the
% functions!

%% First clear/close everything to start fresh: --- USER runs this section ---
clear all; close all; clc

%% setup inputs --- USER edits this section as desired and then runs the section ---
% the main input is a cell array of strings which points to a folder of
% PowerDiva text exports (RLS or DFT) for each subject

saveData=true; % set to true if you want to automatically save the output of rcaSweep into the high-level directory where the data are
printFigures=true; % set to true if you want to automatically print the figures and save them into the high-level directory where the data are
dataLocation='/Volumes/Denali_4D2/rca/'; % change to your data directory
pathnames={[dataLocation,'s001/'],[dataLocation,'s002/'],[dataLocation,'s003/']}; % change folder names to those of your subjects' folders
saveFileName = 'rcaData'; % filename for the .mat file that will be saved in your data directory after analysis finishes

binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse=[1 3 5]; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
condsToUse = [1 2]; % if you want to include all conditions, create a vector here listing all condition numbers
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=3; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports

rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

%% call the function --- USER runs this section (no editing necessary) ---

saveFileNamePath = [dataLocation,'/',saveFileName];
launchAnalysis = false;
if isempty(dir([saveFileNamePath,'.mat']))
    launchAnalysis = true;
else
    fprintf('\n\nWARNING: The file %s.mat already exists for this data directory!\n\n',saveFileName);
    if getYN('Would you like to run rcaSweep anyway and OVERWRITE this file?')
        launchAnalysis = true;
    else
        if getYN('Would you like to run rcaSweep and save the output using a NEW filename?')
            newFileName = input(sprintf('Please type in a new filename and press enter (file already on disk = %s):  ',saveFileName),'s');
            saveFileNamePath = [dataLocation,'/',newFileName];
            launchAnalysis = true;
        end
    end
end
if launchAnalysis
    [rcaData,W,A,covData,noiseData,ozData,ozNoiseData,rcaSettings]=rcaSweep(pathnames,binsToUse,freqsToUse,condsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle);
    if printFigures, print('-dpsc',[saveFileNamePath,'.ps'],'-append'), end
    if saveData, save(saveFileNamePath,'rcaData','W','A','covData','noiseData','ozData','ozNoiseData','rcaSettings'), end
else
    load(saveFileNamePath);
end

%% Examples of provided plotting routines --- USER uncomments/comments plotting calls & runs as desired ---

%% For SNR figures, run ONE of the following:

%[snrFigNums] = plotSnr(rcaData,noiseData,rcaSettings); % generates two plots: one with subplots for each RC, one with subplots for each frequency
%[snrFigNums] = plotSnr(rcaData,noiseData,rcaSettings,[],ozData,ozNoiseData); % additionally creates a plot comparing RC1 to the comparison sensor

% OR provide the name of the comparison electrode as a field in
% plotSettings and then call plotSnr:
plotSettings.comparisonName = 'OZ';
[snrFigNums] = plotSnr(rcaData,noiseData,rcaSettings,plotSettings,ozData,ozNoiseData);

%% For PowerDiva style plots with noise and signal indicated by separate
% markers, with separate subplots for every combination of RC x frequency
plotSettings.conditionColors = [0 0 1; 1 0 0];
plotSettings.showConditions = true;
plotSettings.errorType = 'SEM';
% evaluate ONE of the following:

%figNum = plotFreqByComp(rcaData,noiseData,rcaSettings,plotSettings) % generates one plot without the comparison data
figNum = plotFreqByComp(rcaData,noiseData,rcaSettings,plotSettings,ozData,ozNoiseData); % generates the same plot with an extra column for the comparison sensor

%% Example plot of threshold fit (### want to make fewer requirements for user to set plotSettings..)
plotSettings.titleOn = true;
plotSettings.titleToUse = [];
plotSettings.xlabel = 'Relative Disparity (arc min)';
plotSettings.xTick = rcaSettings.binLevels{1}([1 5 10]);
plotSettings.ymax = [];
plotRcData(rcaData,noiseData,rcaSettings,plotSettings,1,1,1)

%% Print figures to a file
allFigNumbers = [snrFigNums,figNum]; % should include all the returned figure numbers from plotting calls that were already run
if printFigures, for f=allFigNumbers, figure(f); print('-dpsc',[saveFileNamePath,'.ps'],'-append'), end, end
