function [flip_list, corr_list] = componentComparison(axx_in)
    % [flip_list, corr_list] = componentComparison(axx_in)
    
    sensor_axx = cell2mat(permute(axx_in.sensorWave(1,:), [1,3,2]));
    rca_axx = cell2mat(permute(axx_in.rcaWave(1,:), [1,3,2]));
    
    rca_sum = sum(sum(rca_axx, 2), 3);
    sensor_sum = sum(sum(sensor_axx, 2), 3);
    
    rca_corr = corr(rca_sum, sensor_sum);
    
    flip_idx = zeros(1, size(rca_axx,2));
    
    outcomes = logical(dec2bin(0:2^6-1)-'0');
    
    for c = 1:size(outcomes, 1)
        temp_rca = rca_axx;
        temp_rca(:, outcomes(c,:), :) = rca_axx(:,outcomes(c,:),:)*-1;
        temp_sum = sum(sum(temp_rca, 2), 3);
        temp_corr(c) = corr(temp_sum, sensor_sum);
    end
    [corr_list, corr_idx] = sort(temp_corr, 'descend');
    flip_list = outcomes(corr_idx,:);
end

