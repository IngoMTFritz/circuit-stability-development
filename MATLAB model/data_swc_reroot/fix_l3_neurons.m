%% fix disjoint swc files due to axon-dend split in NAVIS

% Use XYZ coordinates of corrupted swc files from meshes
% as input for MST algorithm.
% Use MST trees to generate fixed SWC files to overlay synapses on top at
% NAVIS
clear all;
clc;

%% loadl cells IDs from
idx_fix_swc = [19, 21];

%% load swc files
swcfiles   = dir('*.swc');

opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [8, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "format", "file", "VarName5", "Var6", "Var7"];
opts.SelectedVariableNames = ["format", "file", "VarName5"];
opts.VariableTypes = ["char", "char", "double", "double", "double", "char", "char"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var6", "Var7"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var6", "Var7"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["format", "file", "VarName5"], "ThousandsSeparator", ",");


original_dir = pwd;

for i = idx_fix_swc
    disp(i)
    
    % Import the data
    tree_xyz = readtable(swcfiles(i).name,opts);% Get XYZ coordinates and Volume
    new_tree = MST_tree (1, tree_xyz{:,1}./1000, tree_xyz{:,2}./1000,...% nm to micrometer
        tree_xyz{:,3}./1000, 0, 10, [], [], 'none');
    figure;hold on
    plot_tree     (new_tree, [], [], [], [], '-3l');
    scatter3(new_tree.X(1),new_tree.Y(1),new_tree.Z(1));
    original_dir = pwd;
    
    % Move one directory up
    cd('..');
    % Return to the original directory
    new_tree.X(:) = new_tree.X(:).*1000;
    new_tree.Y(:) = new_tree.Y(:).*1000;
    new_tree.Z(:) = new_tree.Z(:).*1000;

    
    swc_tree (new_tree,append('data_swc_processed/',swcfiles(i).name))
    cd(original_dir);
end

