function cellData=exportToRcaReady(dataPath,how, removeEyes,nanArtifacts,censorEvents)

%This function reorganize the data for 1 subject. By default, we analyze by
%subjects. CellData will be 1Xncondition cell array, within each cell, you
%will have timesample X channels X ntrials. If analyzing by scenes,
%cellData will be a nScenes X nconditions cell array for each subject. The
%cellData will be concatenated later for every subject. 


%% include dot-update filter removal

 
if nargin<5, censorEvents = []; end;
if nargin<4, nanArtifacts = 1; end;
if nargin<3, removeEyes = 1; end;
if nargin<2 || isempty(how),how.nScenes =1;end;

nScenes = how.nScenes;

if how.nScenes >1
    
    nCnd = size(how.allCnd,1);
    %Number of conditions will decide the number of columns in variable "instanceCounter"
    instanceCounter = [(1:nScenes)' zeros(nScenes,nCnd)];

    %instanceCounter count how many times a scene has appeared in each
    %condition. The counted index will be used as the (:,:,n) n here in the
    %third dimension (which is the dimension for trails), so that the trials
    %where the same scene has appeared will be concatenated together.
    
end

RTsegFiles = dir2(fullfile(dataPath, 'RTSeg_*.mat'));
segCount = 0;



for z = 1:numel(RTsegFiles)
    RTfile = fullfile(dataPath, RTsegFiles(z).name);
    load(RTfile); % load data
    if size(TimeLine, 1) > 0;
        segCount = segCount + 1;
        
        
        if isempty(censorEvents) % censor events if requested
            TimeLine = TimeLine(setdiff(1:size(TimeLine, 1), censorEvents));
        else
        end
 
        
         %% read in data trial by trial
        
        for t = 1:size(TimeLine, 1) 

            rawFileName = ['Raw_c' num2str(TimeLine(t).cndNmb, '%03.0f') '_t' num2str(TimeLine(t).trlNmb, '%03.0f')];
            load(fullfile(dataPath, rawFileName));
            rawtrial = double(RawTrial).*repmat(Ampl.', size(RawTrial, 1), 1) + repmat(Shift.' ,size(RawTrial, 1), 1);
            
            if nanArtifacts %make epochs with Artifacts to nan
                    [rowE, colE] = find(IsEpochOK == 0);
                    for idx = 1:length(rowE)
                        rowIdx = ((rowE(idx) - 1)*round(FreqHz) + 1):((rowE(idx) - 1)*round(FreqHz) + round(FreqHz));
                        rawtrial(rowIdx,colE(idx)) = NaN;
                    end
             end
            
            if nScenes == 1 %if analyze by subjects, rawdata will be a 1 X condition cell array, within each cell is the
                            % timeSample X channels X trials 3D matrix
                
                try
                    numNaNs(TimeLine(t).cndNmb,TimeLine(t).trlNmb) = length(find(IsEpochOK == 0));
                    rawdata{1,TimeLine(t). cndNmb}(:, :, TimeLine(t).trlNmb) = rawtrial;
                catch err
                    display(['oops, something went wrong at t = ' num2str(t) ' cnd= ' num2str(TimeLine(t).cndNmb) ' trial= ' num2str(TimeLine(t).trlNmb)]);
                end          
                
            else %if analyze by Scenes
                
                
                try
                    numNaNs(TimeLine(t).cndNmb,TimeLine(t).trlNmb) = length(find(IsEpochOK == 0));
                    whichInstance = find(instanceCounter(:,1)==TimeLine(t).instanceNmb); %find out which instance is in this trial
                    instanceCounter(whichInstance,TimeLine(t).cndNmb+1) = instanceCounter(whichInstance,TimeLine(t).cndNmb+1)+1; %update the count for this instance
                    nInstance = instanceCounter(whichInstance,TimeLine(t).cndNmb+1);%nInstance will be the index for the third dimension
                    rawdata{TimeLine(t).instanceNmb,TimeLine(t). cndNmb}(:, :, nInstance) = rawtrial;
                    %rawdata will be an instance X condition cell array,
                    %where within each cell, there is the 3D matrix
                    %timeSample X channels X numberOfrepitationOfThisInstance
            
                catch err
                    display(['oops, something went wrong at t = ' num2str(t) ' cnd= ' num2str(TimeLine(t).cndNmb) ' trial= ' num2str(TimeLine(t).trlNmb)]);
                end
                
                
                
            end
            
            
            
            
        end
    end
    
end

%     if length(unique(nCond)) > 1
%         warning('Different runs have different number of conditions!')
%     else
%         nCond = length(unique([condInds{:}]));
%     end

nCnd = size(rawdata, 2);


cellData = cell(nScenes,nCnd);

   

for c = 1:nCnd
    for s = 1:nScenes

        if removeEyes
            nTrialsThisCond = size(rawdata{s,c}, 3);
            for tr = 1:nTrialsThisCond
                dataIn = rawdata{s,c}(:, :, tr);
                X = dataIn(:, 1:end - 2); % data channels
                V = dataIn(:, end - 1:end);  % HEOG/VEOG
                if ~isempty(find(isnan(V), 1))
                    disp('ERROR: Eye channels contain NaNs')
                    return;
                end
                dataOut = (eye(size(X, 1)) - V*pinv(V))*X; %A=pinv(V)*X; % transfer function from eyes to data electrodes
                % If a channel has NaN in any epoch, this step will
                % assign NaNs to every epoch in that block for that channel.
                cellData{s,c}(:, :, tr) = dataOut;
            end
        else
            cellData{s,c} = rawdata{s,c}(:, 1:end - 2, :);  % just take out EOG reference electrodes from electrode
        end
    end
end

%     % try to remove prelude and postlude
%     if sum(diff([CndTiming(1).preludeDurSec,CndTiming(1).postludeDurSec,CndTiming(1).stepDurSec]))~=0
%         error(0,'Prelude, Postlude and Bin duration are different!')
%     else
%         binLength = round(CndTiming(1).stepDurSec*FreqHz);
%     end
%     dataIdx = (binLength+1):(size(dataOut,1)-binLength); % assumes that prelude are the same across trials and conditons
%     dataTrim = dataOut(dataIdx,:,:,:);



%     nBins = size(dataTrim,1)/binLength;
%     interData = reshape(dataTrim,binLength,nBins,size(dataTrim,2),size(dataTrim,3),size(dataTrim,4));
%     permData = permute(interData,[1,3,5,2,4]);
%     readyData = permData(:,:,:,:);
%     cellReady = cell(1,size(readyData,3));
%     for r = 1:size(readyData,3)
%         cellReady{r} = squeeze(readyData(:,:,r,:));
%     end
%     clear readyData;
%     readyData = cellReady;
%cd(curDir);
end

