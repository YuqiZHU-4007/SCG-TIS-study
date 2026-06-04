%% Network modularity analysis extracted from original enhanced FC script
% This standalone script keeps the original modularity method unchanged:
%   1. load pre/post matrices from .mat files
%   2. apply Fisher Z transform
%   3. set stats_params.threshold = 1 before modularity analysis
%   4. use the original thresholdMatrix + calculateModularity functions
%   5. compute Pre/Post Q values, relative changes, and statistics

clear; clc; close all;

addpath('C:\安装包\BCT\');
set(0, 'defaultfigurecolor', 'w');

%% ==================== Configuration ====================
labels.path = 'C:\视频任务数据\Network-index\Fuctional connection-Pstatistic\';
labels.operations = {'Stim.10Hz', 'SHAM'};
labels.operations_forload = {'sad_10hz_network.mat', 'sad_sham_network.mat'};
labels.timepoints = {'Pre', 'Post'};

labels.frequencies = { ...
    'Delta (0.5-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-13 Hz)', ...
    'LowBeta (13-20 Hz)', 'HighBeta (20-30 Hz)', 'Gamma1 (30-40 Hz)', ...
    'Gamma2 (40-49 Hz)', 'Gamma3 (51-60 Hz)', 'Gamma4 (60-70 Hz)', ...
    'Gamma5 (70-80 Hz)', 'Gamma6 (80-90 Hz)'};

labels.channels = { ...
    'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6', ...
    'M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6', ...
    'P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};

stats_params.method = 'signrank';
stats_params.correction = 'none';
stats_params.alpha = 0.05;
stats_params.tail = 'auto';
stats_params.savepath = 'C:\视频任务数据\Network-index\Fuctional connection-Pstatistic\';
stats_params.threshold = 1;

fprintf('Starting original-method modularity analysis...\n');

%% ==================== Load and preprocess data ====================
labels = addRegionLabels(labels);
raw_data = loadFunctionalConnectivityData(labels);
data = preprocessData(raw_data);
[Z, ~, X, Y, N] = getDataDimensions(data);

%% ==================== Modularity analysis ====================
[modularity_pre, modularity_post, modularity_changes] = analyzeNetworkModularity( ...
    data, labels, stats_params, Z, X, Y, N);

stats_tables = summarizeModularityStatistics( ...
    data, labels, modularity_pre, modularity_post, modularity_changes, stats_params, Z, X);

%% ==================== Save outputs ====================
if ~exist(stats_params.savepath, 'dir')
    mkdir(stats_params.savepath);
end

save(fullfile(stats_params.savepath, 'modularity_original_method.mat'), ...
    'data', 'labels', 'stats_params', ...
    'modularity_pre', 'modularity_post', 'modularity_changes', 'stats_tables');

writetable(stats_tables.within_pre_post, ...
    fullfile(stats_params.savepath, 'modularity_within_pre_post_original_method.csv'));

if ~isempty(stats_tables.between_relative_change)
    writetable(stats_tables.between_relative_change, ...
        fullfile(stats_params.savepath, 'modularity_between_relative_change_original_method.csv'));
end

fprintf('\nOriginal-method modularity analysis complete.\n');
fprintf('Results saved to: %s\n', stats_params.savepath);

%% ==================== Local functions ====================
function labels = addRegionLabels(labels)
labels.region_mapping = {
    'Fp1', 'Frontal';
    'Fpz', 'Frontal';
    'Fp2', 'Frontal';
    'F7', 'Frontal';
    'F3', 'Frontal';
    'Fz', 'Frontal';
    'F4', 'Frontal';
    'F8', 'Frontal';
    'FC5', 'Fronto-Central';
    'FC1', 'Fronto-Central';
    'FC2', 'Fronto-Central';
    'FC6', 'Fronto-Central';
    'M1', 'Temporal';
    'T7', 'Temporal';
    'C3', 'Central';
    'Cz', 'Central';
    'C4', 'Central';
    'T8', 'Temporal';
    'M2', 'Temporal';
    'CP5', 'Centro-Parietal';
    'CP1', 'Centro-Parietal';
    'CP2', 'Centro-Parietal';
    'CP6', 'Centro-Parietal';
    'P7', 'Parietal';
    'P3', 'Parietal';
    'Pz', 'Parietal';
    'P4', 'Parietal';
    'P8', 'Parietal';
    'POz', 'Parieto-Occipital';
    'O1', 'Occipital';
    'Oz', 'Occipital';
    'O2', 'Occipital'};

labels.regions = unique(labels.region_mapping(:,2))';
labels.channel_regions = cell(size(labels.channels));

for i = 1:length(labels.channels)
    idx = find(strcmp(labels.region_mapping(:,1), labels.channels{i}));
    if ~isempty(idx)
        labels.channel_regions{i} = labels.region_mapping{idx,2};
    else
        labels.channel_regions{i} = 'Unknown';
    end
end
end

function data = loadFunctionalConnectivityData(labels)
Z = min(2, length(labels.operations));
T = length(labels.timepoints);
X = min(11, length(labels.frequencies));

fprintf('Loading data:\n');
fprintf('  operations: %d\n', Z);
fprintf('  timepoints: %d\n', T);
fprintf('  frequencies: %d\n', X);

data = struct();

for z = 1:Z
    file_path = fullfile(labels.path, labels.operations_forload{z});
    if ~exist(file_path, 'file')
        error('Data file not found: %s', file_path);
    end

    a = load(file_path);
    vars = fieldnames(a);

    if any(strcmp(vars, 'pre')) && any(strcmp(vars, 'post'))
        pre = getfield(a, 'pre');
        post = getfield(a, 'post');
        pre_post = cat(5, pre, post);

        Y_z = size(pre_post, 1);
        N = min(32, size(pre_post, 3));

        fprintf('  %s: Y=%d, N=%d\n', labels.operations{z}, Y_z, N);

        data(z).operation_name = labels.operations{z};
        data(z).Y = Y_z;
        data(z).timepoints = struct();

        for t = 1:T
            data(z).timepoints(t).timepoint_name = labels.timepoints{t};
            data(z).timepoints(t).freq = struct();

            for x = 1:X
                data(z).timepoints(t).freq(x).frequency_name = labels.frequencies{x};
                matrices = squeeze(pre_post(1:Y_z, x, 1:N, 1:N, t));
                data(z).timepoints(t).freq(x).matrices = matrices;
            end
        end
    else
        error('Data file %s does not contain pre and post variables.', labels.operations_forload{z});
    end
end
end

function preprocessed_data = preprocessData(data)
[Z, T, X, Y] = getDataDimensions(data);
preprocessed_data = data;

fprintf('Preprocessing: Fisher Z transform...\n');
for z = 1:Z
    for t = 1:T
        for x = 1:X
            Y_z = Y(z);
            for y = 1:Y_z
                matrix = squeeze(data(z).timepoints(t).freq(x).matrices(y,:,:));
                matrix = min(max(matrix, -0.999), 0.999);
                preprocessed_data(z).timepoints(t).freq(x).matrices(y,:,:) = atanh(matrix);
            end
        end
    end
end
end

function [Z, T, X, Y, N] = getDataDimensions(data)
Z = length(data);
T = length(data(1).timepoints);
X = length(data(1).timepoints(1).freq);

Y = zeros(1, Z);
for z = 1:Z
    if isfield(data(z), 'Y')
        Y(z) = data(z).Y;
    else
        if ndims(data(z).timepoints(1).freq(1).matrices) == 3
            Y(z) = size(data(z).timepoints(1).freq(1).matrices, 1);
        else
            Y(z) = 1;
        end
    end
end

if ndims(data(1).timepoints(1).freq(1).matrices) == 3
    N = size(data(1).timepoints(1).freq(1).matrices, 2);
else
    N = size(data(1).timepoints(1).freq(1).matrices, 1);
end

fprintf('Data dimensions:\n');
for z = 1:Z
    fprintf('  operation %d (%s): T=%d, X=%d, Y=%d, N=%d\n', ...
        z, data(z).operation_name, T, X, Y(z), N);
end
end

function [modularity_pre, modularity_post, modularity_changes] = analyzeNetworkModularity(data, labels, stats_params, Z, X, Y, N)
fprintf('Analyzing network modularity with original method...\n');
fprintf('  threshold = %g\n', stats_params.threshold);

modularity_pre = cell(Z, X);
modularity_post = cell(Z, X);
modularity_changes = struct();

for z = 1:Z
    fprintf('  operation %d/%d: %s, subjects=%d\n', z, Z, data(z).operation_name, Y(z));

    for x = 1:X
        pre_vals = zeros(Y(z), 1);
        post_vals = zeros(Y(z), 1);

        for y = 1:Y(z)
            pre_matrix = squeeze(data(z).timepoints(1).freq(x).matrices(y,:,:));
            post_matrix = squeeze(data(z).timepoints(2).freq(x).matrices(y,:,:));

            try
                threshold = stats_params.threshold;
                [pre_mod, ~] = calculateModularity(pre_matrix, threshold, labels);
                [post_mod, ~] = calculateModularity(post_matrix, threshold, labels);

                pre_vals(y) = pre_mod;
                post_vals(y) = post_mod;
            catch ME
                pre_vals(y) = NaN;
                post_vals(y) = NaN;
                fprintf('    %s subject %d failed: %s\n', labels.frequencies{x}, y, ME.message);
            end
        end

        modularity_pre{z, x} = pre_vals;
        modularity_post{z, x} = post_vals;

        valid_idx = ~isnan(pre_vals) & ~isnan(post_vals);
        pre_vals_valid = pre_vals(valid_idx);
        post_vals_valid = post_vals(valid_idx);

        if ~isempty(pre_vals_valid) && ~isempty(post_vals_valid)
            pre_vals_adj = pre_vals_valid;
            pre_vals_adj(abs(pre_vals_valid) < 1e-10) = 1e-10 * sign(pre_vals_valid(abs(pre_vals_valid) < 1e-10));
            pre_vals_adj(pre_vals_adj == 0) = 1e-10;

            modularity_changes(z, x).relative_changes = ...
                (post_vals_valid - pre_vals_adj) ./ abs(pre_vals_adj) * 100;
            modularity_changes(z, x).mean_change = mean(modularity_changes(z, x).relative_changes, 'omitnan');
            modularity_changes(z, x).std_change = std(modularity_changes(z, x).relative_changes, 'omitnan');
            modularity_changes(z, x).pre_mean = mean(pre_vals_valid, 'omitnan');
            modularity_changes(z, x).post_mean = mean(post_vals_valid, 'omitnan');
        else
            modularity_changes(z, x).relative_changes = [];
            modularity_changes(z, x).mean_change = 0;
            modularity_changes(z, x).std_change = 0;
            modularity_changes(z, x).pre_mean = 0;
            modularity_changes(z, x).post_mean = 0;
        end
    end
end

fprintf('Network modularity analysis complete.\n');
end

function [Q, communities] = calculateModularity(connectivity_matrix, threshold, labels)
% Original calculateModularity function. The labels input is kept for API
% compatibility with the original script.

connectivity_matrix_thresh = thresholdMatrix(connectivity_matrix, threshold);
connectivity_matrix_thresh = connectivity_matrix_thresh - diag(diag(connectivity_matrix_thresh));
connectivity_matrix_thresh = max(connectivity_matrix_thresh, connectivity_matrix_thresh');

if all(connectivity_matrix_thresh(:) <= 0)
    Q = 0;
    communities = ones(size(connectivity_matrix_thresh, 1), 1);
    return;
end

M = 1:size(connectivity_matrix_thresh, 1);
Q0 = -1;
Q = 0;
K = 1;
communities = ones(size(connectivity_matrix_thresh, 1), 1);

while Q - Q0 > 1e-5
    K = K + 1;
    Q0 = Q;
    for kk = 1:100
        [community, Q1] = community_louvain(connectivity_matrix_thresh, [], M);
        if Q1 > Q
            communities = community;
        end
        Q = max(Q1, Q);
    end
end
end

function thresholded_matrix = thresholdMatrix(matrix, threshold_percent)
% Original thresholdMatrix function.

triu_indices = triu(true(size(matrix)), 1);
values = matrix(triu_indices);

num_values = length(values);
threshold_idx = round((1 - threshold_percent) * num_values);

if threshold_idx > 0
    sorted_values = sort(values, 'descend');
    threshold_value = sorted_values(threshold_idx);
else
    threshold_value = min(values);
end

thresholded_matrix = matrix;
thresholded_matrix(matrix < threshold_value) = 0;
thresholded_matrix = max(thresholded_matrix, thresholded_matrix');
end

function stats_tables = summarizeModularityStatistics(data, labels, modularity_pre, modularity_post, modularity_changes, stats_params, Z, X)
fprintf('\n=== Modularity statistics ===\n');

operation = {};
frequency = {};
preQ = [];
postQ = [];
relativeChangeMean = [];
relativeChangeStd = [];
pValue = [];
significant = [];

for z = 1:Z
    fprintf('\nWithin-operation Pre vs Post: %s\n', data(z).operation_name);
    fprintf('%-25s %-12s %-12s %-12s %-10s\n', 'Frequency', 'Pre Q', 'Post Q', 'Change %', 'p');
    fprintf('%s\n', repmat('-', 75, 1));

    for x = 1:X
        pre_data = modularity_pre{z, x};
        post_data = modularity_post{z, x};
        pre_data = pre_data(~isnan(pre_data));
        post_data = post_data(~isnan(post_data));

        if ~isempty(pre_data) && ~isempty(post_data)
            test_params = stats_params;
            [p_val, h] = performStatisticalTest(pre_data, post_data, test_params);
        else
            p_val = 1;
            h = 0;
        end

        fprintf('%-25s %-12.4f %-12.4f %-12.2f %-10.4f\n', ...
            labels.frequencies{x}, ...
            modularity_changes(z, x).pre_mean, ...
            modularity_changes(z, x).post_mean, ...
            modularity_changes(z, x).mean_change, ...
            p_val);

        operation{end+1, 1} = data(z).operation_name;
        frequency{end+1, 1} = labels.frequencies{x};
        preQ(end+1, 1) = modularity_changes(z, x).pre_mean;
        postQ(end+1, 1) = modularity_changes(z, x).post_mean;
        relativeChangeMean(end+1, 1) = modularity_changes(z, x).mean_change;
        relativeChangeStd(end+1, 1) = modularity_changes(z, x).std_change;
        pValue(end+1, 1) = p_val;
        significant(end+1, 1) = logical(h);
    end
end

within_pre_post = table(operation, frequency, preQ, postQ, ...
    relativeChangeMean, relativeChangeStd, pValue, significant);

between_relative_change = table();

if Z >= 2
    frequency = {};
    stimChangeMean = [];
    shamChangeMean = [];
    diffShamMinusStim = [];
    pValue = [];
    significant = [];

    fprintf('\nBetween-operation relative change: Stim.10Hz vs SHAM\n');
    fprintf('%-25s %-15s %-15s %-15s %-10s\n', ...
        'Frequency', 'Stim change %', 'SHAM change %', 'SHAM-Stim', 'p');
    fprintf('%s\n', repmat('-', 90, 1));

    for x = 1:X
        stim_changes = modularity_changes(1, x).relative_changes;
        sham_changes = modularity_changes(2, x).relative_changes;
        stim_changes = stim_changes(~isnan(stim_changes));
        sham_changes = sham_changes(~isnan(sham_changes));

        if ~isempty(stim_changes) && ~isempty(sham_changes)
            test_params = stats_params;
            test_params.method = 'ranksum';
            [p_val, h] = performStatisticalTest(stim_changes, sham_changes, test_params);
        else
            p_val = 1;
            h = 0;
        end

        stim_mean = mean(stim_changes, 'omitnan');
        sham_mean = mean(sham_changes, 'omitnan');
        diff_val = sham_mean - stim_mean;

        fprintf('%-25s %-15.2f %-15.2f %-15.2f %-10.4f\n', ...
            labels.frequencies{x}, stim_mean, sham_mean, diff_val, p_val);

        frequency{end+1, 1} = labels.frequencies{x};
        stimChangeMean(end+1, 1) = stim_mean;
        shamChangeMean(end+1, 1) = sham_mean;
        diffShamMinusStim(end+1, 1) = diff_val;
        pValue(end+1, 1) = p_val;
        significant(end+1, 1) = logical(h);
    end

    between_relative_change = table(frequency, stimChangeMean, shamChangeMean, ...
        diffShamMinusStim, pValue, significant);
end

stats_tables.within_pre_post = within_pre_post;
stats_tables.between_relative_change = between_relative_change;
end

function [p_values, h, h_increase, h_decrease, mean_diff] = performStatisticalTest(data1, data2, stats_params)
% Original statistical-test helper, kept to match the original script.

if nargin < 3
    stats_params = struct();
end
if ~isfield(stats_params, 'method'), stats_params.method = 'ranksum'; end
if ~isfield(stats_params, 'tail'), stats_params.tail = 'both'; end

data1 = data1(:);
data2 = data2(:);
data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

if isempty(data1) || isempty(data2)
    p_values = 1;
    h = 0;
    h_increase = 0;
    h_decrease = 0;
    mean_diff = 0;
    return;
end

mean_diff = mean(data2) - mean(data1);

switch lower(stats_params.tail)
    case 'left'
        tail = 'left';
    case 'right'
        tail = 'right';
    case 'auto'
        if mean_diff < 0
            tail = 'right';
        else
            tail = 'left';
        end
    otherwise
        tail = 'both';
end

switch lower(stats_params.method)
    case 'ttest'
        [~, p_values] = ttest(data1, data2, 'Tail', tail);
    case 'ttest2'
        [~, p_values] = ttest2(data1, data2, 'Tail', tail);
    case 'ranksum'
        p_values = ranksum(data1, data2, 'Tail', tail);
    case 'signrank'
        p_values = signrank(data1, data2, 'Tail', tail);
    otherwise
        error('Unsupported statistical method: %s', stats_params.method);
end

if isfield(stats_params, 'alpha')
    alpha = stats_params.alpha;
else
    alpha = 0.05;
end

if strcmpi(stats_params.tail, 'both')
    h = p_values < alpha;
    h_increase = h & (mean_diff > 0);
    h_decrease = h & (mean_diff < 0);
else
    h = p_values < alpha;
    if strcmpi(stats_params.tail, 'left')
        h_increase = 0;
        h_decrease = h;
    else
        h_increase = h;
        h_decrease = 0;
    end
end
end
