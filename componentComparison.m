function [flip_list, corr_list] = componentComparison(input, cond_idx, harm_idx)
    % [flip_list, corr_list] = componentComparison(input)
    
    if isfield(input, 'sensorWave')
        axx_input = true;
        n_cond = size(input.sensorWave, 1);
    else
        axx_input = false;
        n_cond = size(input(1).rca_data, 1);
        n_harm = length(input);
    end
    
    if nargin < 3 || isempty(harm_idx)
        if ~axx_input
            harm_idx = 1:n_harm;
        else
        end
    else
    end
    
    if nargin < 2 || isempty(cond_idx)
        cond_idx = 1:n_cond;
    else
    end
    
    sensor_data = [];
    rca_data = [];
    
    if axx_input
        % concatenate conditions
        for c = 1:length(cond_idx)
            sensor_data = cat(1, sensor_data, cell2mat(permute(input.sensorWave(cond_idx(c),:), [1,3,2])));
            rca_data = cat(1, rca_data,  cell2mat(permute(input.rcaWave(cond_idx(c),:), [1,3,2])));
        end
    else
        comp_idx = contains(input(1).settings.compLabels, 'rc');
        for c = 1:length(cond_idx)
            for h = 1:length(harm_idx)
                sensor_data = cat(1, sensor_data, cell2mat(permute(input(harm_idx(h)).input_data(cond_idx(c),:), [1,3,2])));
                rca_data = cat(1, rca_data,  cell2mat(permute(input(harm_idx(h)).rca_data(cond_idx(c),:), [1,3,2])));
            end
        end
        
    end
    % don't include comparison channel 
    % in rca_data to use for correlation
    rca_data = rca_data(:,comp_idx,:);
    
    % sum over subjects
    rca_sum = nansum(nansum(rca_data, 2), 3);
    sensor_sum = nansum(nansum(sensor_data, 2), 3);
    
    rca_corr = corr(rca_sum, sensor_sum);
    
    flip_idx = zeros(1, size(rca_data,2));
    
    outcomes = logical(dec2bin(0:2^6-1)-'0');
    
    for c = 1:size(outcomes, 1)
        temp_rca = rca_data;
        temp_rca(:, outcomes(c,:), :) = rca_data(:,outcomes(c,:),:)*-1;
        temp_sum = nansum(nansum(temp_rca, 2), 3);
        temp_corr(c) = corr(temp_sum, sensor_sum);
    end
    [corr_list, corr_idx] = sort(temp_corr, 'descend');
    flip_list = outcomes(corr_idx,:);
end

