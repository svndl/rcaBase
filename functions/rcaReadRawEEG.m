function rcaDataOut = rcaReadRawEEG(database)
% Processes subjects rawEEG. Generates and rca-ready matfile   

    rca_path = rca_setPath;
    
    eegRCA = fullfile(rca_path.rcaEEG, database);
    eegSrc = fullfile(rca_path.srcEEG, database);

    proj_dir = eegSrc;

%    proj_dir = uigetdir(eegSrc, 'Select the EEG raw data folder');
%    if (~(proj_dir))
%        error('No directory selected, quitting...');
%    end

    %% take powerdiva export and convert it to cell array format
    list_subj = list_folder(proj_dir);
    
    
    nsubj = numel(list_subj);
    rcaReadyData = fullfile(eegRCA, strcat(database, '_rcaReadyEEG.mat'));
    % check if there is a matfile
    if ~exist(rcaReadyData, 'file')
        
        for s = 1:nsubj
            if (list_subj(s).isdir)
                subjDir = fullfile(proj_dir,  list_subj(s).name);
                try
                    display(['Loading   ' list_subj(s).name]);                  
                    
                    subjEEG = exportToRcaReady(subjDir);
                    rcaData(s, :) = subjEEG(:)';
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
