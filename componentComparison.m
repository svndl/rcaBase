function [flip_list, corr_list] = componentComparison(input, top_sensors, cond_idx, harm_idx, bin_idx)
    % [flip_list, corr_list] = componentComparison(input)
    
    if isfield(input, 'sensorWave')
        axx_input = true;
        n_cond = size(input(1).sensorWave, 1);
        n_sensors = size(input(1).sensorWave{1}, 2);
    else
        axx_input = false;
        n_cond = size(input(1).rca_data, 1);
        n_harm = length(input);
        n_sensors = size(input(1).input_data{1}, 2);
    end
    
    if nargin < 5 || isempty(bin_idx)
        if ~axx_input
            bin_idx = input(1).settings.binIndices;
        else
        end
    else
        if ~axx_input
            if ~ismember(bin_idx, input(1).settings.binIndices)
                msg = 'bin idx must be among the input bins';
                error(msg);
            else
            end
        else
        end
    end
    
    if nargin < 4 || isempty(harm_idx)
        if ~axx_input
            harm_idx = 1:n_harm;
        else
        end
    else
    end
    
    if nargin < 3 || isempty(cond_idx)
        cond_idx = 1:n_cond;
    else
    end
    
    if nargin < 2
        top_sensors = n_sensors;
    else
        if top_sensors > n_sensors
            msg = 'top_sensors cannot be bigger than number of sensors';
            error(msg);
        elseif top_sensors < 1
            msg = 'top_sensors must be at least 1';
            error(msg);
        else
        end
    end
    
    sensor_data = [];
    rca_data = [];
    
    if axx_input
        % concatenate conditions
        for c = 1:length(cond_idx)
            sensor_data = cat(1, sensor_data, cell2mat(permute(input.sensorWave(cond_idx(c),:), [1,3,2])));
            rca_data = cat(1, rca_data,  cell2mat(permute(input.rcaWave(cond_idx(c),:), [1,3,2])));
        end
        % take mean across subjects, then root mean squared
        sensor_rms = rms(nanmean(sensor_data, 3),1);
        [~, sort_idx] = sort(sum(sensor_rms,1), 'descend');
        disp('hello');
    else
        bin_idx = repmat(bin_idx == input(1).settings.binIndices, 2, 1);
        comp_idx = contains(input(1).settings.compLabels, 'rc');
        for h = 1:length(harm_idx)
            % average over trials 
            input(harm_idx(h)).input_data = cellfun(@(x) nanmean(x,3), input(harm_idx(h)).input_data, 'uni', false);
            input(harm_idx(h)).rca_data = cellfun(@(x) nanmean(x,3), input(harm_idx(h)).rca_data, 'uni', false);
            for c = 1:length(cond_idx)
                cur_sensor = cell2mat(permute(input(harm_idx(h)).input_data(cond_idx(c),:), [1,3,2]));
                cur_sensor = cur_sensor(bin_idx, :, :);
                sensor_data = cat(1, sensor_data, cur_sensor);
                cur_rca = cell2mat(permute(input(harm_idx(h)).rca_data(cond_idx(c),:), [1,3,2]));
                cur_rca = cur_rca(bin_idx, :, :);
                rca_data = cat(1, rca_data,  cur_rca);
            end
        end
        % don't include comparison channel 
        % in rca_data to use for correlation
        rca_data = rca_data(:,comp_idx,:);
        
        % select most responsive electrodes
        n_idx = length(harm_idx)*length(cond_idx);
        real_mu = nanmean(sensor_data(1:n_idx,:, :),3);
        imag_mu = nanmean(sensor_data(n_idx+1:end,:, :),3);
        vector_mu = sqrt(real_mu.^2+imag_mu.^2);
        [~, sort_idx] = sort(sum(vector_mu,1), 'descend');
    end
    % grab the most responsive electrodes
    sensor_data = sensor_data(:, sort_idx(1:top_sensors), :);

    % sum over subjects
    rca_sum = nansum(rca_data, 2);
    sensor_sum = nansum(sensor_data, 2);
    
    flip_idx = zeros(1, size(rca_data,2));
    
    n_comp = size(rca_data, 2);
    outcomes = logical(dec2bin(0:2^n_comp-1)-'0');
    
    for c = 1:size(outcomes, 1)
        temp_rca = rca_data;
        temp_rca(:, outcomes(c,:), :) = rca_data(:,outcomes(c,:),:)*-1;
        temp_sum = nansum(temp_rca, 2);
        temp_corr(c) = corr(temp_sum(:), sensor_sum(:));
    end
    [corr_list, corr_idx] = sort(temp_corr, 'descend');
    flip_list = outcomes(corr_idx,:);
end

