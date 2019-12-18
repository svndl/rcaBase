function axxStrct = rcaWaveProject(path_names, W, conds_to_use, type)
% axxStrct = rcaWaveProject(path_names, [W, conds_to_use, type])
% 
% Make axxRCA struct that includes the Axx Wave data 
% and the W (ch x component) matrix from fullRCA
%
% INPUTS:
%   pathnames (required): cell vector of string directory names housing Axx_c00x.mat
%                         or conds x subjects cell matrix of sensor waveforms
%                         (like that returned by this fxn in axxStrct.sensorWave) 
%
%   conds_to_use: vector of conditions to use [ all conditions in data set ]
%   type: prefix 'x' on Axx export folder 'x_Exp_MATL_HCN_128_Avg' [no prefix ]
%   sensorWave: sensorWave data 
% 
% OUTPUTS:
%   struct with the following fields
%       -sensorWave: Wave data from Sensors (loaded from axx-files) 
%       -rcaWave:    Wave data projected through W

    if nargin < 3
        conds_to_use = [];
    else
    end
   
    if nargin < 4
        type = '';
    else
        if ~contains(type, '_')
            type = sprintf('%s_', type);
        else
        end
    end
    
    if any(ischar(path_names{1}))
        path_empty = cell2mat(cellfun(@(x) isempty(exist(x, 'dir')), path_names, 'uni', false));
        if any(path_empty)
            msg = sprintf('\n path %s does not exist', path_names(path_empty));
            error(msg);
        else
        end
        sensorWave = [];
    else
        fprint('\n non-string input given, assuming this is sensor data');
        input_size = size(path_names);
        if numel(input_size) ~= 2
            msg = sprintf('\n sensor input must be two-dimensional \n');
            error(msg);
        else
        end
        if ~iscell(path_names)
            msg = sprintf('\n sensor input must be cell matrix \n');
            error(msg);
        else
        end
        sensorWave = path_names;
    end

    if isempty(sensorWave)
        for f = 1:length(path_names)
            if isempty(conds_to_use)
                axx_matfiles = subfiles(sprintf('%s/%sExp_MATL_HCN_128_Avg/Axx_c*.mat',fileparts(path_names{f}), type), 1);
                axx_matfiles = axx_matfiles(cell2mat(cellfun(@(x) ~contains(x, 'trials'), axx_matfiles, 'uni', false)));
                conds_to_use = cellfun(@(x) str2num(x(end-4)), axx_matfiles);
            else
            end
            axx_matfiles = arrayfun(@(x) ... 
                    sprintf('%s/%sExp_MATL_HCN_128_Avg/Axx_c%03d.mat',fileparts(path_names{f}),type,x), conds_to_use,'uni',false);
            tmpStrct = cellfun(@(x) load(x),axx_matfiles,'uni',false);
            readyStrct = cellfun(@(x) mrC.axx(x),tmpStrct,'uni',false);
            axxStrct.sensorWave(:,f) = cellfun(@(x) x.Wave, readyStrct,'uni',false);
        end
    else
       if isempty(conds_to_use)
            axxStrct.sensorWave = sensorWave;
       else
           axxStrct.sensorWave = sensorWave(conds_to_use, :);
       end
       clear sensorWave;
    end
    
    axxStrct.rcaWave = cell(size(axxStrct.sensorWave));
    
    for z = 1:size(W, 2)
        axxStrct.rcaWave = cellfun(@(x, y) cat(2, x, y), ... 
            axxStrct.rcaWave, rcaProject(axxStrct.sensorWave, W(:,z)), 'uni', false);
    end
end


