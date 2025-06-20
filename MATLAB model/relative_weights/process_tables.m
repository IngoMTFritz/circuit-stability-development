% Load the table from the CSV file
table_data = readtable('l1_adjacency_matrices_norm.csv', 'ReadRowNames', true');

% Get the size of the table
[nRows, nCols] = size(table_data);

% Convert the table to an array for manipulation
data_array = table2array(table_data);

% Iterate over each odd row
for i = 1:2:nRows-1
    % Process odd columns
    for j = 1:2:nCols % For odd columns
        if j < nCols % Ensure j+1 exists
            mean_value = mean([data_array(i, j), data_array(i+1, j+1)]);
            data_array(i, j) = mean_value;       % Replace (i, j)
            data_array(i+1, j+1) = mean_value;  % Replace (i+1, j+1)
        end
    end
    % Process even columns
    for k = 2:2:nCols % For even columns
        if k > 1 % Ensure k-1 exists
            mean_value = mean([data_array(i, k), data_array(i+1, k-1)]);
            data_array(i, k) = mean_value;       % Replace (i, k)
            data_array(i+1, k-1) = mean_value;  % Replace (i+1, k-1)
        end
    end
end

% Convert the modified array back to a table
modified_table = array2table(data_array, 'RowNames', table_data.Properties.RowNames, 'VariableNames', table_data.Properties.VariableNames);

% Save the modified table to a new CSV file
writetable(modified_table, 'l1_adjacency_matrices_norm_mean.csv', 'WriteRowNames', true);


%% L3

% Load the table from the CSV file
table_data = readtable('l3_adjacency_matrices_norm.csv', 'ReadRowNames', true');

% Get the size of the table
[nRows, nCols] = size(table_data);

% Convert the table to an array for manipulation
data_array = table2array(table_data);

% Iterate over each odd row
for i = 1:2:nRows-1
    % Process odd columns
    for j = 1:2:nCols % For odd columns
        if j < nCols % Ensure j+1 exists
            mean_value = mean([data_array(i, j), data_array(i+1, j+1)]);
            data_array(i, j) = mean_value;       % Replace (i, j)
            data_array(i+1, j+1) = mean_value;  % Replace (i+1, j+1)
        end
    end
    % Process even columns
    for k = 2:2:nCols % For even columns
        if k > 1 % Ensure k-1 exists
            mean_value = mean([data_array(i, k), data_array(i+1, k-1)]);
            data_array(i, k) = mean_value;       % Replace (i, k)
            data_array(i+1, k-1) = mean_value;  % Replace (i+1, k-1)
        end
    end
end

% Convert the modified array back to a table
modified_table = array2table(data_array, 'RowNames', table_data.Properties.RowNames, 'VariableNames', table_data.Properties.VariableNames);

% Save the modified table to a new CSV file
writetable(modified_table, 'l3_adjacency_matrices_norm_mean.csv', 'WriteRowNames', true);
