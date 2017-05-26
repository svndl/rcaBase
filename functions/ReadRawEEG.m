function rcaDataOut = ReadRawEEG(database,how)
% Processes subjects rawEEG.    
% if how.nScenes = 1 , the rcaDataOut will be a subject X condition cell array, if
% how.nScenes > 1, the rcaDataOut will be a scene X condition cell array
% default setting is nScenes =1

if nargin<2 || isempty(how), how.nScenes = 1; end

    natSc_path = natSc_setPath(database,how);
    nScenes = how.nScenes;
    eegRCA = fullfile(natSc_path.rcaEEG, database);
    if (~exist(eegRCA, 'dir'));
        mkdir(eegRCA);
    end;
    eegSrc = fullfile(natSc_path.srcEEG, database);

    proj_dir = eegSrc;

%    proj_dir = uigetdir(eegSrc, 'Select the EEG raw data folder');
%    if (~(proj_dir))
%        error('No directory selected, quitting...');
%    end

    %% take powerdiva export and convert it to cell array format
    list_subj = list_folder(proj_dir);
    
    
    nsubj = numel(list_subj);
    if nScenes == 1
        rcaReadyData = fullfile(eegRCA, strcat(database, '_rcaReadyEEG','_bySubjects','.mat'));
    else
        rcaReadyData = fullfile(eegRCA, strcat(database, '_rcaReadyEEG','_byScenes','.mat'));
    end
    % check if there is a matfile
    if ~exist(rcaReadyData, 'file')
        
        for s = 1:nsubj
            if (list_subj(s).isdir)
                subjDir = fullfile(proj_dir,  list_subj(s).name);
                try
                    display(['Loading   ' list_subj(s).name]);                  
                    
                    subjEEG = exportToRcaReady(subjDir,how);
                    
                    
                    if nScenes ==1
                    rcaData(s, :) = subjEEG(:);
                    else
                        if s ==1
                            rcaData = subjEEG;
                        else
                            %Concatenate the repetition of the same
                            %instance from different subjects
                            rcaData = cellfun(@(a,b) cat(3,a,b), rcaData, subjEEG,'UniformOutput',false);
                        end
                    end
                    
                catch
                    display(['Warning, could not load   ' list_subj(s).name]);                                   
                %do nothing
                end
                
            end
        end
        rcaDataOut = reshape(rcaData(~cellfun('isempty', rcaData)), size(rcaData));
        
        save(rcaReadyData, '-v7.3', 'rcaDataOut');   
    else
        load(rcaReadyData);
    end
end
