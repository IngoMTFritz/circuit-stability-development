% Define bin edges from 0 to 1 in steps of 0.1
function bins = plotBinnedAverages(X_clean, Y_clean, range)
    % Define bin edges
    bin_edges = range;
    
    % Initialize arrays to store the average of each bin
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    bin_averages = zeros(size(bin_centers));
    
    % Loop through each bin and calculate the average Y value for the bin
    for i = 1:length(bin_edges)-1
        % Find indices of X_clean that fall within the current bin
        bin_indices = X_clean >= bin_edges(i) & X_clean < bin_edges(i+1);
        
        % Calculate the average Y value for the current bin
        if any(bin_indices)
            bin_averages(i) = mean(Y_clean(bin_indices));
        else
            bin_averages(i) = NaN; % Handle empty bins
        end
    end
    
    % Plot the average Y values for each bin
    plot(bin_centers, bin_averages, '--','LineWidth', 2);
    bins = bin_averages ;

end

