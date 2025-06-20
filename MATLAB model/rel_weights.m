function modified_values = rel_weights(ti, matrix_file)
    % Function to process a neuron name, find the corresponding column in a matrix,
    % and return the modified values (multiplied by 100 and ceiled).
    %
    % Arguments:
    %   ti: A structure with a .name field containing the neuron name.
    %   matrix_file: Path to the matrix CSV file.
    %
    % Returns:
    %   modified_values: A column vector of modified values.

    % Extract the neuron name and remove "_dendrite"
    neuron_name = ti.name;
    cleaned_name = erase(neuron_name, '_dendrite');
    cleaned_name = erase(cleaned_name, 'l1_');
    cleaned_name = erase(cleaned_name, 'l3_');
    
    % Load the matrix from the specified file
    matrix_table = readtable(matrix_file, 'ReadRowNames', true);
    headers = readcell(matrix_file, 'Range', '1:1'); % Read first row as headers
    headers  = headers(2:end);
    % Find the corresponding column
    col_idx = find(strcmp(headers, cleaned_name), 1);
    
    if isempty(col_idx)
        error('Neuron name "%s" not found in the matrix headers.', cleaned_name);
    end

    % Retrieve the column values
    column_values = matrix_table(:, col_idx);
    column_values = column_values{:,:};

    % Multiply by 100 and apply ceil
    modified_values = column_values * 100;
end
