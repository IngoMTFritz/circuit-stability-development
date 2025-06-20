%% Use real synapse locations and density for distributed injections
%  Use different average radius for l1 and l3
%  Gm from Gynay et al and the same in l1 and l3
%
clear;
clc;

c_l1 = [0.59, 0.78, 0.87];
c_l3 = [0.39, 0.35, 0.64];
ddaC = [0.87, 0.32, 0.22];
vada = [1.00, 0.58, 0.25];
vdaB = [1.00, 0.84, 0.48];

%% convert to micrometers, give average diameter and passive properties from Gunay plos CB

% Specify the folder containing your SWC files
folderName = 'data_swc_processed';

% Get a list of all SWC files in the folder
fileList = dir(fullfile(folderName, '*.swc'));

% Preallocate a cell array to hold the loaded trees
t = cell(length(fileList), 1);%length(fileList)

% Loop through each file and load it using load_tree
mD = 0.78; % mean diameter for aCC neuron cell in Gunay

% load syn densities
data = readmatrix('data/density_l1_l3_dendrite.csv');
last_column = data(:, end);

% Compute the mean and standard deviation
mean_den= mean(last_column);
std_den = std(last_column);
min_den = min(last_column);
max_den = max(last_column);

%L1
for i = 1:12
    disp(i)
    filePath = fullfile(folderName, fileList(i).name);
    tree = load_tree(filePath);
    tree.X(:) = tree.X(:)/1000;
    tree.Y(:) = tree.Y(:)/1000;
    tree.Z(:) = tree.Z(:)/1000;
    tree.D(:) = mD; % mean diameter for l1 KC
    tree.Gm = 3.796e-5; % generic values from Gunay
    tree.Ri = 41.47; % generic values from Gunay
    tree = resample_tree(tree,0.1); % average density in l1 and l3
    t{i} = tree;
end

% L3 different radius than l1?
for i = 13:length(t)
    disp(i)
    filePath = fullfile(folderName, fileList(i).name);
    tree = load_tree(filePath);
    tree.X(:) = tree.X(:)/1000;
    tree.Y(:) = tree.Y(:)/1000;
    tree.Z(:) = tree.Z(:)/1000;
    tree.D(:) = mD; % mean diameter for l1 KC
    tree.Gm = 3.796e-5; % generic values from Gunay
    tree.Ri = 41.47; % generic values from Gunay
    tree = resample_tree(tree,0.1); % average density in l1 and l3
    t{i} = tree;
end

save('data/t', 't');
%% Plot l1

Vsyn            = 10*1/mean_den; %mV generic spiking threshold in Gunay

%store values
synprox_mean_store = zeros(length(t), 50); % 50 relative weights activated
synprox_std_store  = zeros(length(t), 50);
syndist_mean_store = zeros(length(t), 50);
syndist_std_store  = zeros(length(t), 50);
nodes_all_store    = zeros(length(t), 50);


for i = 1:length(t)
    disp(i)
    tree             = t{i};
    N                = length (tree.X);
    Pvec             = Pvec_tree (tree);
    % Given or known constants
    Ee               = 40;% generic fly cell from Cuntz et al.
    Ei               = 0;
    Gm               = tree.Gm;      % S/cm2
    GmU              = 10 * Gm;      % in nS/um2
    Idist            = Vsyn*t{1}.Gm * pi * mD*10;                            % pA/(0.01*1/mean_denum)
    gsynT            = Idist / Ee;
    ge_total         = gsynT / 1000;                          % Total ge when homogeneous
    Vdist            = (Idist/(1/mean_den) /10) ./ (t{1}.Gm * pi * mD+ (gsynT/(1/mean_den)/10));
    gi               = zeros  (N, 1);
    I                = zeros  (N, 1);
    L    (i)         = sum  (len_tree (t{i})) ;    % in um
    n                = floor  (linspace (0, round(mean_den * L(i)), 50));
    nn               = length (n);
    RR               = zeros  (nn, 1);
    synprox          = zeros  (nn, 1);
    syndist          = zeros  (nn, 1);
    synact           = zeros  (nn, 1); % active synapses
    
    for counter = 1:nn
        syns = round(mean_den * L(i));  % Number of synaptic entries to modify
        
        syndist = zeros(100,1); % Store 100 trials per counter
        synprox = zeros(100,1); % Store 100 trials per counter
        synact  = zeros(100,1); % Store node count per trial
            
        for trial = 1:100
            
            % Generate log-normal distributed ge values
            ge_ln = lognrnd(0, 1, syns, 1);
            ge_ln = ge_ln / sum(ge_ln) * (syns * ge_total);
            
            % Initialize ge with zeros
            ge = zeros(N, 1);
            
            % Select 'syns' linearly spaced indices and update them
            lin_indices = round(linspace(1, N, syns)); % Get 'syns' evenly spaced indices
            ge(lin_indices) = ge_ln; % Add (ge_ln) only to these entries
            
            R = n(counter); % Set threshold R
            
            % ------- Processing gedist -------
            gedist = ge;
            [~, P] = sort(Pvec); % Sort Pvec in ascending order
            
            % Find the first R elements within lin_indices
            valid_indices = intersect(P, lin_indices, 'stable'); % Keep order from P
            if length(valid_indices) >= R
                gedist(valid_indices(1:R)) = 0;  % Zero only first R valid indices
            end
            
            syn = syn_tree(tree, gedist, gi, Ee, Ei, I, 'none');
            syndist(trial) = syn(1);
            
            % ------- Processing geprox -------
            geprox = ge;
            [~, P] = sort(Pvec, 'descend'); % Sort Pvec in descending order
            
            % Find the first R elements within lin_indices
            valid_indices = intersect(P, lin_indices, 'stable'); % Keep order from P
            if length(valid_indices) >= R
                geprox(valid_indices(1:R)) = 0;  % Zero only first R valid indices
            end
            
            syn = syn_tree(tree, geprox, gi, Ee, Ei, I, 'none');
            synprox(trial) = syn(1);
            
            % Count active synapses
            synact(trial) = length(find(geprox));
        end
        
        % Store the statistics for each counter
            synprox_mean_store(i, counter) = mean(synprox);
            synprox_std_store(i, counter)  = std(synprox);
            syndist_mean_store(i, counter) = mean(syndist);
            syndist_std_store(i, counter)  = std(syndist);            
            nodes_all_store(i, counter)    = mean(synact);
    end

end


%% Plot final voltage

c                    = [0.4940 0.1840 0.5560];
hc                   = [0.25 0.3  0.25];
figure;
clf;
hold             on;


scatter(L(1:12), synprox_mean_store(1:12,1), 10, c_l1, 'filled');
scatter(L(13:end), synprox_mean_store(13:end,1), 10, c_l3, 'filled');
scatter(L(1:12), syndist_mean_store(1:12,1), 10, c_l1, 'filled');
scatter(L(13:end), syndist_mean_store(13:end,1), 10, c_l3, 'filled');

pagesize         = [5 4]; % (x, y) size in cm

yline(Vdist, '--','linewidth',0.5); % Add label and line style

xlabel  ('Length (\mum)');
ylabel  ('Voltage (mV)');
ylim([6 10]);
title ('Mean Model')
set(gcf,'renderer','Painters')

set              (gca, ...
    'ytick',                   0 : 2 : 10, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'xscale',                  'log', ...
    'fontname',                'arial');

nRMSE   = sqrt (mean (((synprox_mean_store(:,1) - Vdist)./synprox_mean_store(:,1)).^2)) % 0.12% 

tprint           (          ...
    'output_plots/Figure 5A - L vs V',   ...
    '-HR -eps -jpg',        ...
    pagesize);



%% Plot Rel weigths l1

syn_mean = mean([synprox_mean_store(1:12,:);syndist_mean_store(1:12,:)], 1);
syn_std = sqrt( mean([synprox_std_store(1:12,:).^2 + synprox_mean_store(1:12,:).^2; ...
                       syndist_std_store(1:12,:).^2 + syndist_mean_store(1:12,:).^2], 1) ...
               - syn_mean.^2);

figure;
hold             on;
psyn = flip(100 * n./R); % Relative weights [%]

errorbar(psyn, syn_mean, syn_std, 'LineWidth', 1,...
    'DisplayName', 'Distal Removal','CapSize', 0, 'Color', c_l1);

% Target Vdist as reference
xax = 0:100;
yax = Vdist * xax / 100;
plot(xax, yax, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Target V_{dist}');

xlim             ([0 100]);
ylim             ([0 10]);
%pagesize         = [2 2]; % (x, y) size in cm
pagesize         = [5 4]; % (x, y) size in cm

xlabel  ('Relative Weights [%]');
ylabel  ('Voltage [mV]');
title   ('L1');

set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 2 : 10, ...
    'xtick',                   0 : 25 : 100, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

rw = repmat(psyn, size(synprox_mean_store(1:12,:),1), 1); % relative weights
compute_r2([rw(:);rw(:)],[synprox_mean_store(1:12,:),synprox_mean_store(1:12,:)],'n')


tprint           (          ...
    'output_plots/Figure 5B - Rel Weights L1',   ...
    '-HR -eps -jpg',        ...
    pagesize);


%% Plot abs weights L1

figure;
hold             on;

l1_abs = nodes_all_store(1:12,:);
l1_prox= synprox_mean_store(1:12,:);
l1_dist= syndist_mean_store(1:12,:);

scatter(l1_abs(:), l1_prox(:), 6, c_l1, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(l1_abs(:), l1_dist(:), 6, c_l1, 'filled', 'MarkerFaceAlpha', 0.7);

ylim([0 10])

xlabel  ('# Synapses');
ylabel  ('Voltage [mV]');
pagesize         = [5 4]; % (x, y) size in cm


set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 2 : 12, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

compute_r2([l1_abs(:);l1_abs(:)],[l1_prox(:);l1_dist(:)],'y')
bin_averages = plotBinnedAverages([l1_abs(:);l1_abs(:)], [l1_prox(:);l1_dist(:)], 0:20:300);

tprint           (          ...
    'output_plots/Figure 5B - Abs Weights L1',   ...
    '-HR -eps -jpg',        ...
    pagesize);


%%%%%%%%%%%%%%%%%%%%%%
%% Plot Rel weigths l3
%%%%%%%%%%%%%%%%%%%%%%

syn_mean = mean([synprox_mean_store(13:end,:);syndist_mean_store(13:end,:)], 1);
syn_std = sqrt( mean([synprox_std_store(13:end,:).^2 + synprox_mean_store(13:end,:).^2; ...
                       syndist_std_store(13:end,:).^2 + syndist_mean_store(13:end,:).^2], 1) ...
               - syn_mean.^2);

figure;
hold             on;
psyn = flip(100 * n./R); % Relative weights [%]

errorbar(psyn, syn_mean, syn_std, 'LineWidth', 1,...
    'DisplayName', 'Distal Removal','CapSize', 0,'Color', c_l3);

% Target Vdist as reference
xax = 0:100;
yax = Vdist * xax / 100;
plot(xax, yax, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Target V_{dist}');

xlim             ([0 100]);
ylim             ([0 10]);
%pagesize         = [2 2]; % (x, y) size in cm
pagesize         = [5 4]; % (x, y) size in cm

xlabel  ('Relative Weights [%]');
ylabel  ('Voltage [mV]');
title   ('L3');

set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 2 : 12, ...
    'xtick',                   0 : 50 : 100, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

rw = repmat(psyn, size(synprox_mean_store(13:end,:),1), 1); % relative weights
compute_r2([rw(:);rw(:)],[synprox_mean_store(13:end,:),syndist_mean_store(13:end,:)],'n')



tprint           (          ...
    'output_plots/Figure 5C - Rel Weights L3',   ...
    '-HR -eps -jpg',        ...
    pagesize);


%% Plot abs weights L3

figure;
hold             on;

l3_abs = nodes_all_store(13:end,:);
l3_prox= synprox_mean_store(13:end,:);
l3_dist= syndist_mean_store(13:end,:);

scatter(l3_abs(:), l3_prox(:), 6, c_l3, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(l3_abs(:), l3_dist(:), 6, c_l3, 'filled', 'MarkerFaceAlpha', 0.7);

ylim([0 10]);
xlim([0 1100]);

xlabel  ('# Synapses');
ylabel  ('Voltage [mV]');
pagesize         = [5 4]; % (x, y) size in cm


set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 2 : 10, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

compute_r2([l3_abs(:);l3_abs(:)],[l3_prox(:);l3_dist(:)],'y')
bin_averages = plotBinnedAverages([l3_abs(:);l3_abs(:)], [l3_prox(:);l3_dist(:)], 0:100:1200);

tprint           (          ...
    'output_plots/Figure 5C - Abs Weights L3',   ...
    '-HR -eps -jpg',        ...
    pagesize);


%%%%%%%%%%%%%%%%%%%%%
%% Plot Abs L1-L3
%%%%%%%%%%%%%%%%%%%%%


figure;
hold             on;

scatter(l3_abs(:), l3_prox(:), 6, c_l3, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(l3_abs(:), l3_dist(:), 6, c_l3, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(l1_abs(:), l1_prox(:), 6, c_l1, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(l1_abs(:), l1_dist(:), 6, c_l1, 'filled', 'MarkerFaceAlpha', 0.7);

ylim([0 10]);
xlim([0 1100]);
pagesize         = [5 4]; % (x, y) size in cm

xlabel  ('# Synapses');
ylabel  ('Voltage [mV]');


set(gcf,'renderer','Painters')
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75], ...
    'ytick',                   0 : 2 : 10, ...
    'XMinorTick',              'on', ...
    'YMinorTick',              'on', ...
    'ticklength',              [0.096 0.24] ./ max (pagesize), ...
    'tickdir',                 'out', ...
    'linewidth',               0.5, ...
    'fontsize',                6, ...
    'fontname',                'arial');

compute_r2([nodes_all_store(:);nodes_all_store(:)],...
                  [synprox_mean_store(:);syndist_mean_store(:)],'y')
bin_averages = plotBinnedAverages([nodes_all_store(:);nodes_all_store(:)],...
                  [synprox_mean_store(:);syndist_mean_store(:)], 0:100:1200);

tprint           (          ...
    'output_plots/Figure 5D - Abs Weights L1 vs L3',   ...
    '-HR -eps -jpg',        ...
    pagesize);


%% Morphology
clf;
hold             on;
filePath = fullfile(folderName, fileList(1).name);
tree = load_tree(filePath);
tree.X(:) = tree.X(:)/1000;
tree.Y(:) = tree.Y(:)/1000;
tree.Z(:) = tree.Z(:)/1000;

hp               = plot_tree     (tree, c_l1, [], [], [], '-3l');

filePath = fullfile(folderName, fileList(13).name);
tree = load_tree(filePath);
tree.X(:) = tree.X(:)/1000;
tree.Y(:) = tree.Y(:)/1000;
tree.Z(:) = tree.Z(:)/1000;
hp               = plot_tree     (tree, c_l3, [], [], [], '-3l');
%set              (hp, ...
%    'edgecolor',               'none');
hp               = line ([0 10], [0 0] + 70);
set              (hp, ...
    'linewidth',               0.5, ...
    'color',                   [0 0 0]);
axis off
pagesize         = [5 5]; % (x, y) size in cm
set              (gca, ...
    'ActivePositionProperty',  'position', ...
    'position',                [0.2 0.2 0.75 0.75]);
tprint           (          ...
    'output_plots/Figure 5B - morphologies',   ...
    '-HR -eps -jpg',        ...
    pagesize);



