%% PAC-Behavior-Mechanism Integrated Analysis Script
clear; clc; close all;
%% ==================== 1. 数据初始化与参数设置 ====================
isstim=1;
savepath='C:\视频任务数据\Network-index\PAC\通道间\';
if isstim==1;
    load('C:\视频任务数据\Network-index\PAC\通道间\10hz_pac_z_m.mat');
elseif isstim==0;
    load('C:\视频任务数据\Network-index\PAC\通道间\sham_pac_z_m.mat');
end
pac_z_score=squeeze(nanmean(pac_z_score,3));
n_sub = size(pac_z_score, 2);
alpha_axis = linspace(0.5,20,39);
gamma_axis = linspace(30,90,60);
load('C:\视频任务数据\videovalence.mat');
videovalence([[1,2,6],20+[1,5]],:) = [];
A = table2array(videovalence);
ii = 10; ind1 = find(A(:,2)==isstim);ind=ind1;
behavior_data = (A(ind,ii)-A(ind,ii-1));%./(A(ind,ii-1))
load('C:\视频任务数据\BDI.mat');
BDI([[1,2,6],20+[1,5]],:) = [];
A = table2array(BDI);
ii = 4; ind1 = find(A(:,2)==isstim);ind=ind1;
behavior_data2 = (A(ind,ii)-A(ind,ii-1));
data_pre=squeeze(pac_z_score(1,:,:,:,:,:));data_post=squeeze(pac_z_score(2,:,:,:,:,:));
all_diff = (data_post - data_pre);
mean_diff_all = squeeze(nanmean(all_diff, 1));
chan_names = {'Fp1','F7','F3','FC5','T7','CP6','O1','Oz','O2'};
n_chan = length(chan_names);
alpha_lvl = 0.05;
smooth_sigma = 1.2;
min_cluster_size=round(0.01*numel(alpha_axis)*numel(gamma_axis));
sig_mask_final = cell(n_chan, n_chan);
h_matrix = zeros(n_chan, n_chan);
row_hub_score = zeros(n_chan, 1);
col_hub_score = zeros(n_chan, 1);
for i = 1:n_chan
    row_p_vals = ones(n_chan, 1);
    for j = 1:n_chan
        curr_pre  = squeeze(data_pre(:, i, j, :, :));
        curr_post = squeeze(data_post(:, i, j, :, :));
        p_map = zeros(size(curr_pre, 2), size(curr_pre, 3));
        for fg = 1:size(curr_pre, 2)
            for fa = 1:size(curr_pre, 3)
                p_map(fg, fa) = signrank(curr_post(:, fg, fa), curr_pre(:, fg, fa));
            end
        end
        binary_mask = p_map < 0.05;corrected_mask = bwareaopen(binary_mask, min_cluster_size);
        row_p_vals(j) = median(p_map(find(corrected_mask==1)));
        sig_mask_final{i, j} = corrected_mask;
    end
    [~, h_row] = fdr(row_p_vals, alpha_lvl);
    h_matrix(i, :) = h_row;
end
%% MI change
chani=[1,3:6];
chanj=7:9;
a_sub = squeeze(nanmean(nanmean(data_pre(:,chanj,chani, :, :), 3), 2));
b_sub = squeeze(nanmean(nanmean(data_post(:,chanj,chani, :, :), 3), 2));
p_avg_map = zeros(size(a_sub, 2), size(a_sub, 3));
for fg = 1:size(a_sub, 2)
    for fa = 1:size(a_sub, 3)
        p_avg_map(fg, fa) = signrank(b_sub(:, fg, fa), a_sub(:, fg, fa));
    end
end
sig_mask_avg = p_avg_map <= 1;
sig_mask_avg = bwareaopen(sig_mask_avg, min_cluster_size);
sig_mask_avg(:,[1:15,26:39])=0;
cluster_idx = find(sig_mask_avg == 1);
n_sub_active = size(a_sub, 1);
sub_cluster_change = zeros(n_sub_active, 1);
for s = 1:n_sub_active
    pre_val  = mean(a_sub(s, cluster_idx));
    post_val = mean(b_sub(s, cluster_idx));
    sub_cluster_change(s) = post_val - pre_val;
end
%% heatmap
chani=[1,3:6];
chanj=7:9;
a_sub = squeeze(nanmean(nanmean(data_pre(:,chani,chanj, :,:), 3), 2));
b_sub = squeeze(nanmean(nanmean(data_post(:,chani,chanj, :, :), 3), 2));
p_avg_map = zeros(size(a_sub, 2), size(a_sub, 3));
for fg = 1:size(a_sub, 2)
    for fa = 1:size(a_sub, 3)
        p_avg_map(fg, fa) = signrank(b_sub(:, fg, fa), a_sub(:, fg, fa));
    end
end
sig_mask_avg = p_avg_map < 0.05;
sig_mask_avg = bwareaopen(sig_mask_avg, min_cluster_size);
for hh=1
    a=squeeze(nanmean(a_sub,1));
    b=squeeze(nanmean(b_sub,1));
    diff_map = b - a; % Post - Pre
    diff_sm = imgaussfilt(diff_map, smooth_sigma);
    h=figure('Color', 'w', 'Position', [10 50 2000 400]);
    data_plots = {a, b};
    titles = {'Pre-Stimulation', 'Post-Stimulation'};
    c_raw = [a(:); b(:)]; clims_raw = [0.02, 0.54];%clims_raw = [min(c_raw), max(c_raw)];
    for k = 1:2
        subplot(1, 5, k);
        imagesc(alpha_axis, gamma_axis, imgaussfilt(data_plots{k}, smooth_sigma), clims_raw);
        hold on; axis xy; colormap(redblue);
        title(titles{k}); xlabel('Phase (Hz)');
        if k == 1, ylabel('Amplitude (Hz)'); end
        colorbar;
    end
    ax3 = subplot(1, 5, 3);
    h_img = imagesc(alpha_axis, gamma_axis, diff_sm);
    hold on; axis xy;
    max_d = max(abs(diff_sm(:)));
    caxis([-max_d, max_d]);
    set(h_img, 'AlphaData', sig_mask_avg * 0.8 + 0.2);
    if any(sig_mask_avg(:))
        contour(alpha_axis, gamma_axis, sig_mask_avg, 1, 'k', 'LineWidth', 1.8);
    end
    [max_v, max_i] = max(diff_map(:) .* sig_mask_avg(:));
    if max_v > 0
        [py, px] = ind2sub(size(diff_map), max_i);
        plot(alpha_axis(px), gamma_axis(py), 'p', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'y', 'MarkerSize', 14);
        text(alpha_axis(px)+0.5, gamma_axis(py), ...
            sprintf('Peak: %.1f/%.1fHz', alpha_axis(px), gamma_axis(py)), ...
            'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.6]);
    end
    title('Statistical Difference (Post - Pre)', 'Color', 'k');
    xlabel('Phase(Hz)');
    cb3 = colorbar; %ylabel(cb3, '\Delta PAC (Z-score)');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    for kk=1:2
        sig_mask_avg(:,[1:15,26:39])=0;
        cluster_idx = find(sig_mask_avg == 1);
        if isempty(cluster_idx)
            warning('未检测到显著簇，请检查统计阈值或簇过滤大小。');continue;
        end
        n_sub_active = size(a_sub, 1);
        sub_cluster_change = zeros(n_sub_active, 1);
        for s = 1:n_sub_active
            pre_val  = mean(a_sub(s, cluster_idx));
            post_val = mean(b_sub(s, cluster_idx));
            sub_cluster_change(s) = post_val - pre_val;
        end
        switch kk
            case 1; behavior_datay=behavior_data;labely='\Delta Arousal';
            case 2; behavior_datay=behavior_data2;labely='\Delta BDI';
        end
        ax4 = subplot(1, 5, kk+3);
        defaults.title = '';
        defaults.xlabel = '\Delta PAC (Post - Pre)';
        defaults.ylabel = labely;
        defaults.show_equation = false;
        [p_val,r_val]=plot_eeg_behavior_correlation(sub_cluster_change, behavior_datay, defaults);
        hold on;axis square;
        text_str = sprintf('r_{s} = %.3f\n p = %.4f', r_val, p_val);
        set(gcf, 'Renderer', 'Painters');
        if p_val < 0.05
            % text(min(sub_cluster_change), max(behavior_datay), ' Significant Association ', ...
            %     'Color', 'w', 'BackgroundColor', [0.4 0.7 0.4], 'FontWeight', 'bold');
        end
    end
end

%% bubblechart of PAC-behav corr
sig_mask_alpha_gamma=[];
sig_mask_alpha_gamma_p={};
sig_mask_alpha_gamma_r={};
for chani=[1,3:6];
    for chanj=7:9;
        a_sub = squeeze(nanmean(nanmean(data_pre(:,chani,chanj, :,:), 3), 2));%[1,3:6],7:9
        b_sub = squeeze(nanmean(nanmean(data_post(:,chani,chanj, :, :), 3), 2));%chani,chanj
        p_avg_map = zeros(size(a_sub, 2), size(a_sub, 3));
        for fg = 1:size(a_sub, 2)
            for fa = 1:size(a_sub, 3)
                p_avg_map(fg, fa) = signrank(b_sub(:, fg, fa), a_sub(:, fg, fa));
            end
        end
        sig_mask_avg = p_avg_map < 0.05;
        sig_mask_avg = bwareaopen(sig_mask_avg, min_cluster_size);
        for hh=1
            a=squeeze(nanmean(a_sub,1));
            b=squeeze(nanmean(b_sub,1));
            diff_map = b - a; % Post - Pre
            diff_sm = imgaussfilt(diff_map, smooth_sigma);
            h=figure('Name', chan_names{chani}, 'Color', 'w', 'Position', [10 50 2000 400]);
            data_plots = {a, b};
            titles = {'Pre-Stimulation', 'Post-Stimulation'};
            c_raw = [a(:); b(:)]; clims_raw = [0.02, 0.54];%clims_raw = [min(c_raw), max(c_raw)];
            for k = 1:2
                subplot(1, 5, k);
                imagesc(alpha_axis, gamma_axis, imgaussfilt(data_plots{k}, smooth_sigma), clims_raw);
                hold on; axis xy; colormap(redblue);
                title(titles{k}); xlabel('Phase (Hz)');
                if k == 1, ylabel('Amplitude (Hz)'); end
                colorbar;
            end
            ax3 = subplot(1, 5, 3);
            h_img = imagesc(alpha_axis, gamma_axis, diff_sm);
            hold on; axis xy;
            max_d = max(abs(diff_sm(:)));
            caxis([-max_d, max_d]);
            set(h_img, 'AlphaData', sig_mask_avg * 0.8 + 0.2);
            if any(sig_mask_avg(:))
                contour(alpha_axis, gamma_axis, sig_mask_avg, 1, 'k', 'LineWidth', 1.8);
            end
            [max_v, max_i] = max(diff_map(:) .* sig_mask_avg(:));
            if max_v > 0
                [py, px] = ind2sub(size(diff_map), max_i);
                plot(alpha_axis(px), gamma_axis(py), 'p', 'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', 'y', 'MarkerSize', 14);
                text(alpha_axis(px)+0.5, gamma_axis(py), ...
                    sprintf('Peak: %.1f/%.1fHz', alpha_axis(px), gamma_axis(py)), ...
                    'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.6]);
            end
            title('Statistical Difference (Post - Pre)', 'Color', 'k');
            xlabel('Phase(Hz)');
            cb3 = colorbar; %ylabel(cb3, '\Delta PAC (Z-score)');
            set(findall(gcf,'-property','FontSize'),'FontSize',14);
            for kk=1:2
                sig_mask_avg(:,[1:15,26:39])=0;
                cluster_idx = find(sig_mask_avg == 1);
                if isempty(cluster_idx)
                    warning('未检测到显著簇，请检查统计阈值或簇过滤大小。');continue;
                end
                n_sub_active = size(a_sub, 1);
                sub_cluster_change = zeros(n_sub_active, 1);
                for s = 1:n_sub_active
                    pre_val  = mean(a_sub(s, cluster_idx));
                    post_val = mean(b_sub(s, cluster_idx));
                    sub_cluster_change(s) = post_val - pre_val;
                end
                switch kk
                    case 1; behavior_datay=behavior_data;labely='\Delta Arousal';
                    case 2; behavior_datay=behavior_data2;labely='\Delta BDI';
                end
                ax4 = subplot(1, 5, kk+3);
                defaults.title = '';
                defaults.xlabel = '\Delta PAC (Post - Pre)';
                defaults.ylabel = labely;
                defaults.show_equation = false;
                [p_val,r_val]=plot_eeg_behavior_correlation(sub_cluster_change, behavior_datay, defaults);
                hold on;axis square;
                text_str = sprintf('r_{s} = %.3f\n p = %.4f', r_val, p_val);
                set(gcf, 'Renderer', 'Painters');
                if p_val < 0.05
                    % text(min(sub_cluster_change), max(behavior_datay), ' Significant Association ', ...
                    %     'Color', 'w', 'BackgroundColor', [0.4 0.7 0.4], 'FontWeight', 'bold');
                end
                sig_mask_alpha_gamma(chani,chanj)=mean(sub_cluster_change);
                sig_mask_alpha_gamma_p{kk}(chani,chanj)=p_val;
                sig_mask_alpha_gamma_r{kk}(chani,chanj)=r_val;
            end
            saveas(h,string(fullfile(savepath,strcat('PAC change_',chan_names(chani),'-',chan_names(chanj),'_with arousal.png'))));
            exportgraphics(h,string(fullfile(savepath,strcat('PAC change_',chan_names(chani),'-',chan_names(chanj),'_with arousal.pdf'))), 'ContentType', 'vector', 'BackgroundColor', 'none');
            close all;
        end
    end
end
chani=[1,3:6];chanj=7:9;
for iii=1:2
    labels1 = chan_names(chani);
    labels2 = chan_names(chanj);
    raw_p_data = sig_mask_alpha_gamma_p{iii}(chani, chanj)';
    raw_p_data (raw_p_data  == 0) = nan;
    data = sig_mask_alpha_gamma(chani, chanj)';
    data2 = sig_mask_alpha_gamma_r{iii}(chani, chanj)';
    [row, col] = ind2sub(size(data), (1:numel(data))');
    h=figure('Color', 'w', 'Units', 'inches', 'Position', [2, 2, 8, 3]); hold on;
    b = bubblechart(col, row, -data(:), data2(:));
    my_map = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)'];
    colormap(redblue);clim([-0.6 0.6]); %flipud(my_map)
    c = colorbar;
    ylabel(c, 'Effect Size', 'FontSize', 10, 'FontWeight', 'bold');
    bubblesize([1 25]);
    blgd = bubblelegend('- R', 'Location', 'northeastoutside', 'Style', 'vertical');
    blgd.NumBubbles = 3;
    blgd.Box = 'off';
    v_p = raw_p_data(:);
    for i = 1:length(v_p)
        if isnan(v_p(i)), continue; end
        str = '';
        if v_p(i) < 0.001
            str = '***';
        elseif v_p(i) < 0.01
            str = '**';
        elseif v_p(i) < 0.05
            str = '*';
        end
        if ~isempty(str)
            text(col(i), row(i), str, ...
                'Color', 'k', ...
                'FontSize', 12, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end
    ax = gca;
    xticks(1:length(labels1)); xticklabels(labels1);
    yticks(1:length(labels2)); yticklabels(labels2);
    ax.YDir = 'reverse';
    ax.XAxisLocation = 'top';
    ax.TickDir = 'out';
    ax.LineWidth = 1.2;
    box off;
    xlim([0.5 length(labels1)+0.5]);
    ylim([0.5 length(labels2)+0.5]);
end
saveas(h,[savepath,'channel PAC and arousal.png']);
exportgraphics(h, [savepath,'channel PAC and arousal.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

%% Topographical distribution of source and target hub contributions
source_idx = [1,3:6]; 
target_idx = 7:9; 
for chani=1
    pathway_row_hub = zeros(n_chan, 1);
    pathway_col_hub = zeros(n_chan, 1);
    for i = source_idx
        for j = target_idx
            if h_matrix(i, j) == 1 && any(any(sig_mask_final{i, j}(:,16:25)))
                area_val = sum(sum(sig_mask_final{i, j}(:,16:25)));
                pathway_row_hub(i) = pathway_row_hub(i) + area_val;
                pathway_col_hub(j) = pathway_col_hub(j) + area_val;
            end
        end
    end
    h=figure('Color', 'w', 'Name', 'Pathway-Specific Hub (Frontal Phase -> Occipital Amp)', 'Position', [100 100 900 450]);
    subplot(1,2,1);
    b1 = bar(source_idx, pathway_row_hub(source_idx), 'FaceColor', [0.2 0.5 0.7]);
    set(gca, 'XTick', source_idx, 'XTickLabel', chan_names(source_idx));
    xtickangle(45);
    title('Phase Hubs (Providers)', 'FontSize', 12);
    ylabel('Significant Area to Occipital (Pixels)');
    grid on; axis square;
    subplot(1,2,2);
    b2 = bar(target_idx, pathway_col_hub(target_idx), 'FaceColor', [0.8 0.4 0.3]);
    set(gca, 'XTick', target_idx, 'XTickLabel', chan_names(target_idx));
    xtickangle(45);
    title('Amplitude Hubs (Receivers)', 'FontSize', 12);
    ylabel('Significant Area from Frontal (Pixels)');
    grid off; axis square;
    [max_row_v, max_row_i] = max(pathway_row_hub(source_idx));
    actual_row_idx = source_idx(max_row_i);
    [max_col_v, max_col_i] = max(pathway_col_hub(target_idx));
    actual_col_idx = target_idx(max_col_i);
end
labels.channels = {'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6',...
    'M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6',...
    'P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};
source_idx=[1,5,9,14,23];contribution1=zeros(1,32);contribution1(source_idx)=pathway_row_hub([1,3:6]);
source_idx=[30:32];contribution2=zeros(1,32);contribution2(source_idx)=pathway_col_hub([7:9]);
chan_locs = create_manual_locs(labels.channels);
for ii=1
    try
        Th = [chan_locs.theta]';
        Rd = [chan_locs.radius]';
        [x_coords, y_coords] = pol2cart((90-Th)*pi/180, Rd);

        if max(Rd) > 0.5
            x_coords = x_coords * (0.5 / max(Rd));
            y_coords = y_coords * (0.5 / max(Rd));
        end

        x_coords = x_coords(:);
        y_coords = y_coords(:);
    catch
        error('坐标计算失败。');
    end
    h=figure('Color', 'w', 'Position', [100, 100, 600, 600]);
    topoplot([], chan_locs, 'style', 'blank', 'electrodes', 'off', 'headrad', 'rim', 'shading', 'interp');
    hold on;
    bubble_sizes = ((contribution1 - min(contribution1(:))) / (max(contribution1(:)) - min(contribution1(:)) + eps)) * 400 + 30;
    scatter(x_coords, y_coords, bubble_sizes(:), contribution1(:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerFaceAlpha', 0.8,'MarkerFaceColor','r');
    bubble_sizes = ((contribution2 - min(contribution2(:))) / (max(contribution2(:)) - min(contribution2(:)) + eps)) * 400 + 30;
    scatter(x_coords, y_coords, bubble_sizes(:), contribution2(:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerFaceAlpha', 0.8,'MarkerFaceColor','b');
    [~, sort_idx] = sort(contribution1, 'descend');
    for i = 1:length(sort_idx)
        curr_idx = sort_idx(i);
        if contribution1(curr_idx) > 0
            text(x_coords(curr_idx), y_coords(curr_idx) + 0.04, ...
                chan_locs(curr_idx).labels, 'FontSize', 11, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'Color', 'b');
        end
    end
    [~, sort_idx] = sort(contribution2, 'descend');
    for i = 1:length(sort_idx)
        curr_idx = sort_idx(i);
        if contribution2(curr_idx) > 0
            text(x_coords(curr_idx), y_coords(curr_idx) + 0.04, ...
                chan_locs(curr_idx).labels, 'FontSize', 11, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'Color', 'b');
        end
    end
    title('Corrected Topo-Bubble Map');
    axis equal; axis tight;
end

saveas(h,[savepath,' Hub 贡献度.png']);
exportgraphics(h, [savepath,'Hub 贡献度.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

