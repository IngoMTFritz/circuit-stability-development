%% Simulations based on relative weights

% the density of synapses is sampled from the distro of observed values

% synaptic strength follows a log normal distro with total current the same
% as if all synapses were homogeneous in strength

% relative weights are normalised and averaged over left and right values

clear;
clc;
load             ('data/t.mat');
c_l1 = [0.59, 0.78, 0.87];
c_l3 = [0.39, 0.35, 0.64];
ddaC = [0.87, 0.32, 0.22];
vada = [1.00, 0.58, 0.25];
vdaB = [1.00, 0.84, 0.48];
%% load syn densities
data = readmatrix('data/density_l1_l3_dendrite.csv');
last_column = data(:, end);

% Compute the mean and standard deviation
mean_den= mean(last_column);
std_den = std(last_column);
min_den = min(last_column);
max_den = max(last_column);


%% Plot relative weights
% Load the tables
tablel1 = readtable('relative_weights/l1_adjacency_matrices_norm_mean.csv', 'ReadRowNames', true);
tablel3 = readtable('relative_weights/l3_adjacency_matrices_norm_mean.csv', 'ReadRowNames', true);

% Convert tables to numeric arrays
vector1 = table2array(tablel1(:,:)); % Extract numerical values from table 1
vector2 = table2array(tablel3(:,:)); % Extract numerical values from table 2

% Reshape into column vectors for scatter plot and correlation computation
vector1 = vector1(:)*100; % Convert to a column vector
vector2 = vector2(:)*100; % Convert to a column vector

% Apply threshold: Set values below 1 to 0 (Tomoko Nature 2015)
%vector1(vector1 < 1) = 0;
%vector2(vector2 < 1) = 0;

% Ensure both vectors have the same length
min_length = min(length(vector1), length(vector2));
vector1 = vector1(1:min_length);
vector2 = vector2(1:min_length);

pagesize         = [5 5]; % (x, y) size in cm

% Create scatter plot
figure; hold on;
s1 = scatter(vector1(1:24), vector2(1:24), 10, ddaC, 'filled'); 
s2 = scatter(vector1(25:48), vector2(25:48), 10, vada, 'filled'); 
s3 = scatter(vector1(49:72), vector2(49:72), 10, vdaB, 'filled'); 

% Add legend with corresponding color names
legend([s1, s2, s3], {'ddaC', 'vada', 'vdaB'}, 'Location', 'best');

xlabel('L1 Relative Weights');
ylabel('L3 Relative Weights');
xlim([-0.5 10])
ylim([-0.5 10])

set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 2 : 12, ...
    'xtick',                   0 : 2 : 12, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

% Compute Pearson correlation coefficient
[r, p] = corr(vector1, vector2, 'Type', 'Pearson');
compute_r2(vector1, vector2,'y')

% Display the correlation results
disp(['Pearson correlation coefficient (r): ', num2str(r)]);
disp(['p-value: ', num2str(p)]);

tprint           (          ...
    'output_plots/Figure 7C - Rel Weights L1vs L3',   ...
    '-HR -eps -jpg',        ...
    pagesize);

%save('data_output/Fig7C_rel_.mat', 'vector1', 'vector2')  % Save variables you want

%% Plot l1

   % Given or known constants
    Ee               = 40;% generic fly cell from Cuntz et al.
    Ei               = 0;
    Gm               = t{1}.Gm;      % S/cm2
    GmU              = 10 * Gm;      % in nS/um2
    delta            = compute_overshoot(mean_den, GmU, mean(t{1}.D), Ee, 8);
    Vdist            = (8+delta)*1/mean_den; % 8mV generic spiking threshold in Gunay

% Vdist estimates the voltage needed at a distal dendritic site to produce an 8 mV local depolarization with synaptic conductances.
% The expression Vdist = (8 + delta=2) * 1 / mean_den approximates the current needed (Idist), assuming I = V * Gmpd.
% The extra 2 mV is added to compensate for the fact that in conductance-based synapses,
% Vsyn = Isyn * Gmpd + Gsyn — the synaptic conductance reduces voltage efficacy, so more current
% (and thus more initial voltage) is required to reach the same 8 mV target.
%
% This formulation also assumes a synaptic density of 1 per unit length. Since our neurons use a lower
% synaptic density (mean_den), each synapse must be stronger to maintain the same local depolarization.
% Therefore, Vdist is scaled inversely with mean_den to ensure that Vsyn remains sufficient (~8 mV).

% Preallocate matrices for storing results (12 mdiv for `i`, 6 columns for realtive weights)
synrand_mean_store = zeros(length(t), 6);
synrand_std_store  = zeros(length(t), 6);
weights_all_store  = zeros(length(t), 6);
nodes_all_store    = zeros(length(t), 6);

for i = 1:length(t)
    disp(i)
    tree             = t{i};
    N                = length(tree.X);
    Pvec             = Pvec_tree(tree);
    mD               = t{1}.D(1);
    
    if i<13
        n  = rel_weights(tree, 'relative_weights/l1_adjacency_matrices_norm_mean.csv'); % relative weights from mdiv->ln L1
    else
        n  = rel_weights(tree, 'relative_weights/l3_adjacency_matrices_norm_mean.csv'); % relative weights from mdiv->ln L1
    end
    
    % Given or known constants
    Idist            = Vdist*GmU * pi * mD;   % pA/(0.01*1/mean_den/um) to get units right
    gsynT            = Idist / Ee; % Ohm's law
    ge_total         = gsynT / 1000; % Total ge when homogeneous. /1000 to match units of syn_tree function
    Vsyn             = (Idist/(1/mean_den)) ./ (GmU * pi * mD+ (gsynT/(1/mean_den)));
    gi               = zeros(N, 1);
    I                = zeros(N, 1);
    L    (i)         = sum  (len_tree (t{i})) ;    % in um

    for den = 1:10 %trials of different densities per neuron
        % Sample mS from a normal distribution
        mS = normrnd(mean_den, std_den);
        if mS < min_den, mS = min_den; end
        if mS > max_den, mS = max_den; end

        syns = round(mS * L(i));  % Number of synaptic entries to modify
        
        for counter = 1:6 % Since `nn = 6`
            syn_trials = zeros(100,1); % Store 10 trials per counter
            node_counts = zeros(100,1); % Store node count per trial
            
            for trial = 1:100
                disp([i, den, counter, trial])
                
                % Generate log-normal distributed ge values
                ge_random = lognrnd(0, 1, syns, 1);
                ge_random = ge_random / sum(ge_random) * (syns * ge_total);
                
                gerand  = zeros(N,1);
                percentage_to_zero = 100 - n(counter);
                num_to_zero = syns-round((percentage_to_zero / 100) * syns);
                [~, iR] = sort(rand(N, 1));
                gerand(iR(1:num_to_zero)) = ge_random(1:num_to_zero); % Set some entries to zero
                
                syn = syn_tree(tree, gerand, gi, Ee, Ei, I, 'none');
                syn_trials(trial) = syn(1); % Store for this trial
                node_counts(trial) = length(find(gerand));
            end
            
            % Store the statistics for each counter
            synrand_mean_store(i, counter) = mean(syn_trials);
            synrand_std_store(i, counter)  = std(syn_trials);
            weights_all_store(i, counter)  = n(counter);
            nodes_all_store(i, counter)    = mean(node_counts);
        end
    end
end


%% Plot Voltages from Relative Weights

% Assume synrand_mean_store and synrand_std_store are both 24x6 matrices

% Extract first and second halves
mean_l1 = synrand_mean_store(1:12, :);
mean_l3 = synrand_mean_store(13:24, :);
std_l1 = synrand_std_store(1:12, :);
std_l3 = synrand_std_store(13:24, :);

% Flatten matrices into vectors for plotting
x = mean_l1(:);  % X values (First Half)
y = mean_l3(:);  % Y values (Second Half)
x_err = std_l1(:);  % X-axis error (First Half STD)
y_err = std_l3(:); % Y-axis error (Second Half STD)

% Create figure
figure; hold on;
xlim([-0.1 1.3]);
ylim([-0.1 1.3]);

% Grouped data (reshaping each subset into column vectors)
x1 = mean_l1(:,1:2); y1 = mean_l3(:,1:2);
x2 = mean_l1(:,3:4); y2 = mean_l3(:,3:4);
x3 = mean_l1(:,5:6); y3 = mean_l3(:,5:6);

% Convert matrices to column vectors
x1 = x1(:); y1 = y1(:);
x2 = x2(:); y2 = y2(:);
x3 = x3(:); y3 = y3(:);

% Standard deviations (reshape into column vectors)
std_x1 = std_l1(:,1:2); std_y1 = std_l3(:,1:2);
std_x2 = std_l1(:,3:4); std_y2 = std_l3(:,3:4);
std_x3 = std_l1(:,5:6); std_y3 = std_l3(:,5:6);

std_x1 = std_x1(:); std_y1 = std_y1(:);
std_x2 = std_x2(:); std_y2 = std_y2(:);
std_x3 = std_x3(:); std_y3 = std_y3(:);

% Scatter plot for different groups
s1 = scatter(x1, y1, 10, ddaC, 'filled'); % ddaC color
s2 = scatter(x2, y2, 10, vada, 'filled'); % vada color
s3 = scatter(x3, y3, 10, vdaB, 'filled'); % vdaB color

%save('data_output/Fig7D_vol_ddaC.mat', 'x1', 'y1');  % Save variables you want
%save('data_output/Fig7D_vol_vada.mat', 'x2', 'y2');  % Save variables you want
%save('data_output/Fig7D_vol_vbaB.mat', 'x3', 'y3');  % Save variables you want

% Add error bars for each group
% errorbar(x1, y1, std_y1, std_y1, std_x1, std_x1, 'o', ...
%     'MarkerSize', 2, 'MarkerEdgeColor', ddaC, 'MarkerFaceColor', ddaC, 'CapSize', 0, 'Color', ddaC);
% 
% errorbar(x2, y2, std_y2, std_y2, std_x2, std_x2, 'o', ...
%     'MarkerSize', 2, 'MarkerEdgeColor', vada, 'MarkerFaceColor', vada, 'CapSize', 0, 'Color', vada);
% 
% errorbar(x3, y3, std_y3, std_y3, std_x3, std_x3, 'o', ...
%     'MarkerSize', 2, 'MarkerEdgeColor', vdaB, 'MarkerFaceColor', vdaB, 'CapSize', 0, 'Color', vdaB);

% Add legend
legend([s1, s2, s3], {'ddaC', 'vada', 'vdaB'}, 'Location', 'best');

% Labels and Title
xlabel('Voltage L1 (mV)');
ylabel('Voltage L3 (mV)');

pagesize         = [5 5]; % (x, y) size in cm

set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 0.4 : 2, ...
    'xtick',                   0 : 0.4 : 2, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');
    
[r, p] = corr(mean_l1(:), mean_l3(:), 'Type', 'Pearson');
disp(r)
compute_r2(mean_l1, mean_l3,'y')


tprint           (          ...
    'output_plots/Figure 7D - Voltage L1 vs L3',   ...
    '-HR -eps -jpg',        ...
    pagesize);

%% Plot L1 and L3 voltage responses vs relative weights

figure; hold on;
rel_l1 = weights_all_store(1:12, :);
rel_l3 = weights_all_store(13:24, :);
abs_l1 = nodes_all_store(1:12, :);
abs_l3 = nodes_all_store(13:24, :);


figure; hold on;
errorbar(rel_l1, mean_l1, std_l1, 'Color',c_l1 ,'MarkerSize', 2.5,...
         'DisplayName', 'L1', 'CapSize', 0, 'LineStyle', 'none', 'Marker', 'o');
errorbar(rel_l3, mean_l3, std_l3, 'Color',c_l3 ,'MarkerSize', 2.5, ...
         'DisplayName', 'L3', 'CapSize', 0, 'LineStyle', 'none', 'Marker', 'o');
     
%save('data_output/Fig7E_rel_L1.mat', 'rel_l1', 'std_l1', 'mean_l1');  % Save variables you want
%save('data_output/Fig7E_rel_L3.mat', 'rel_l3', 'std_l3', 'mean_l3');  % Save variables you want


xlim             ([-0.1 10.5]);
ylim             ([-0.1 1.6]);
%pagesize         = [2 2]; % (x, y) size in cm
pagesize         = [5 4]; % (x, y) size in cm

xlabel  ('Relative Weights [%]');
ylabel  ('Voltage [mV]');

set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 0.2 : 4, ...
    'xtick',                   0 : 2 : 12, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

[r, p] = corr(weights_all_store(:),synrand_mean_store(:), 'Type', 'Pearson');
compute_r2(weights_all_store(:),synrand_mean_store(:),'y')


tprint           (          ...
    'output_plots/Figure 7E - L1_L3 Voltage vs Rel Weights',   ...
    '-HR -eps -jpg',        ...
    pagesize);


%% Plot L1 and L3 voltage responses vs num of synapses

figure; hold on;
errorbar(abs_l1, mean_l1, std_l1, 'Color',c_l1 ,'MarkerSize', 2.5,...
         'DisplayName', 'L3', 'CapSize', 0, 'LineStyle', 'none', 'Marker', 'o');
errorbar(abs_l3, mean_l3, std_l3, 'Color',c_l3 ,'MarkerSize', 2.5,...
         'DisplayName', 'L1', 'CapSize', 0, 'LineStyle', 'none', 'Marker', 'o');
     
%save('data_output/Fig7E_abs_L1.mat', 'abs_l1', 'std_l1', 'mean_l1');  % Save variables you want
%save('data_output/Fig7E_abs_L3.mat', 'abs_l3', 'std_l3', 'mean_l3');  % Save variables you want
     

xlim             ([-0.1 80]);
ylim             ([-0.1 1.6]);
pagesize         = [5 4]; % (x, y) size in cm

xlabel  ('#Syns');
ylabel  ('Voltage [mV]');

set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 0.2 : 4, ...
    'xtick',                   0 : 10 : 100, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

[r, p] = corr(nodes_all_store(:),synrand_mean_store(:), 'Type', 'Pearson');
compute_r2(nodes_all_store(13:24,:),synrand_mean_store(13:24,:),'y')


tprint           (          ...
    'output_plots/Figure 7E - L1_L3 Voltage - Abs Weights',   ...
    '-HR -eps -jpg',        ...
    pagesize);

