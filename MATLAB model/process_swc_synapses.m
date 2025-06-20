% Define the input and output directories
input_folder = 'data';
output_folder = 'data_swc_processed';

% Create the output folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Get a list of all SWC files in the input folder
swc_files = dir(fullfile(input_folder, '*.swc'));

% Loop over each SWC file
for i = 1:length(swc_files)
    % Get the full path of the current SWC file
    input_file = fullfile(input_folder, swc_files(i).name);
    
    % Read the content of the SWC file
    file_content = fileread(input_file);
    
    % Split the file into lines
    lines = strsplit(file_content, '\n');
    
    % Initialize a cell array to hold the processed lines
    processed_lines = {};
    
    % Process each line
    for j = 1:length(lines)
        line = strtrim(lines{j});
        
        % Skip empty lines
        if isempty(line)
            continue;
        end
        
        % Check if the line is a comment
        if startsWith(line, '#')
            processed_lines{end+1} = line; 
            continue;
        end
        
        % Process data lines
        tokens = strsplit(line);
        if numel(tokens) >= 8
            % Replace column 2 with column 8, then remove column 8
            tokens{2} = tokens{8};
            tokens(8) = [];
            % Convert the updated tokens back to a string
            processed_line = strjoin(tokens, ' ');
            processed_lines{end+1} = strtrim(processed_line); %#ok<AGROW>
        else
            % Keep the line as it is if it doesn't have enough columns
            processed_lines{end+1} = line; %#ok<AGROW>
        end
    end
    
    % Join the processed lines into a single string
    processed_content = strjoin(processed_lines, '\n');
    
    % Write the processed content to a new file in the output folder
    output_file = fullfile(output_folder, swc_files(i).name);
    fid = fopen(output_file, 'w');
    if fid == -1
        error('Cannot write to file: %s', output_file);
    end
    fwrite(fid, processed_content);
    fclose(fid);
end

disp('Processing complete.');
