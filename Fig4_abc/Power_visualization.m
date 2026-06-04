function main()
% Pre-Post,Pos-Neg,Individual,Frqe,factor,channel
% clc;clear all;
% loadpath='C:\视频任务数据\论文数据整理\task2 data and code\Fig4_abc\Results\';
% load([loadpath,'\sham_relative_psd_power.mat']);
% psd_data=squeeze(nanmean(relative_power(:,2,:,:,:,:),5));size(psd_data);
% save([loadpath,'\sham_psd_neg.mat']);
% load([loadpath,'\10hz_relative_psd_power.mat'])
% psd_data=squeeze(nanmean(relative_power(:,2,:,:,:,:),5));size(psd_data);
% save([loadpath,'\10hz_psd_neg.mat']);
Analysis_Script();
end

function Analysis_Script()
clear; clc; close all;
set(0, 'DefaultfigureColor', 'w');
labels = generateLabels();
chanlocs = create_manual_locs(labels.channels);
stim_path = fullfile(labels.path, labels.operations_forload{1});
sham_path = fullfile(labels.path, labels.operations_forload{2});
load(stim_path); Stim_Raw = psd_data;
load(sham_path); Sham_Raw = psd_data;
load('C:\视频任务数据\videovalence.mat');
videovalence([[1,2,6],20+[1,5]],:) = [];
A = table2array(videovalence);
ii = 10; ind1 = find(A(:,2)==1);ind0 = find(A(:,2)==0);
behavior_stim1 = (A(ind1,ii)-A(ind1,ii-1));
behavior_sham1 = (A(ind0,ii)-A(ind0,ii-1));
load('C:\视频任务数据\BDI.mat');
BDI([[1,2,6],20+[1,5]],:) = [];
A = table2array(BDI);
ii = 4; ind1 = find(A(:,2)==1);ind0 = find(A(:,2)==0);
behavior_stim2 = (A(ind1,ii)-A(ind1,ii-1));A(ind1,:);
behavior_sham2 = (A(ind0,ii)-A(ind0,ii-1));A(ind0,:);
[~, n_subs, n_freqs, n_chans] = size(Stim_Raw);
Stim_Delta = squeeze((Stim_Raw(2,:,:,:) - Stim_Raw(1,:,:,:))./ abs(Stim_Raw(1,:,:,:)));
Sham_Delta =squeeze((Sham_Raw(2,:,:,:) - Sham_Raw(1,:,:,:))./ abs(Sham_Raw(1,:,:,:)));

%% === topoplot and  Bar plot ===
sig_mask=[];isboxplot=1;
yll=[[-20,20];[-10 10];[-8 8];[-8 8];[-10 10];[-15 10];[-15 10];[-15 15];[-15 10];[-15 15];[-15 15]];
for target_band_idx = 1:n_freqs
    target_band_name = labels.frequencies{target_band_idx};
    fprintf('\n正在分析特定频段: %s\n', target_band_name);
    Stim_Delta_Band = squeeze(Stim_Delta(:, target_band_idx, :));
    Sham_Delta_Band = squeeze(Sham_Delta(:, target_band_idx, :));
    t_values = zeros(1, n_chans);
    p_values = zeros(1, n_chans);
    for ch = 1:n_chans
        [p_inc,~ , stats_inc] = ranksum(Stim_Delta_Band(:, ch), Sham_Delta_Band(:, ch),'tail','right');
        [p_dec,~ , stats_dec] = ranksum(Stim_Delta_Band(:, ch), Sham_Delta_Band(:, ch),'tail','left');
        p_values(ch) = min(p_inc,p_dec); 
        z_inc = 0; z_dec = 0;
        if isfield(stats_inc, 'zval'), z_inc = stats_inc.zval; end
        if isfield(stats_dec, 'zval'), z_dec = stats_dec.zval; end
        t_values(ch) = max(z_inc, z_dec);
    end
    [p_values,~] = fdr(p_values,0.05);
    sig_mask(:,:,target_band_idx) = p_values < 0.05;
    if true %target_band_idx == 2 || target_band_idx == 3
        h=figure('Color', 'w', 'Name', ['TI Effect Topo - ' target_band_name], 'Position', [100, 100, 1600, 600]);
        subplot(1, 3, 1);
        topoplot(mean(Stim_Delta_Band, 1), chanlocs, 'style', 'map', 'electrodes', 'off', 'headrad', 'rim', 'shading', 'interp');
        colormap('jet');hc=colorbar;setcolorbar(hc);%'maplimits', yll(target_band_idx,:),
        title({'Stim Group', target_band_name}, 'FontSize', 12); set(gca,'FontName','Arial');
        subplot(1, 3, 2);
        topoplot(mean(Sham_Delta_Band, 1), chanlocs, 'style', 'map', 'electrodes', 'off', 'headrad', 'rim', 'shading', 'interp');
        colormap('jet');colormap('jet');hc=colorbar;setcolorbar(hc);
        title({'Sham Group', target_band_name}, 'FontSize', 12); set(gca,'FontName','Arial');
        subplot(1, 3, 3); yl(1)=round(min(t_values)); yl(2)=round(max(t_values));
        topoplot(t_values, chanlocs, 'maplimits',[yl], 'style', 'map','electrodes', 'off', 'headrad', 'rim', 'shading', 'interp',...
            'emarker2', {find(sig_mask(:,:,target_band_idx)==1), 'o', 'w', 6, 1});hold on;
        colormap('jet');colormap('jet');hc=colorbar;setcolorbar(hc);
        title({'Interaction (Z-values)', target_band_name}, 'FontSize', 12, 'FontWeight', 'bold');set(gca,'FontName','Arial');
        saveas(h,[labels.path,['TI Effect Topo - ' target_band_name],'.png']);
        exportgraphics(h, [labels.path,['TI Effect Topo - ' target_band_name],'.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
        close(h);
        h=figure('Color', 'w', 'Name', ['All Channels Detail - ' target_band_name], 'Position', [50, 50, 1800, 1000]);
        grid_rows = 4; grid_cols = 8;
        for ch = 1:n_chans
            subplot(grid_rows, grid_cols, ch);
            data_stim = Stim_Delta_Band(:, ch);
            data_sham = Sham_Delta_Band(:, ch);
            m_stim = mean(data_stim); sem_stim = std(data_stim)/sqrt(length(data_stim));
            m_sham = mean(data_sham); sem_sham = std(data_sham)/sqrt(length(data_sham));
            if isboxplot
                color_stim = [255,128,128]/255;
                color_sham = [141,141,255]/255;
                group_names = {char(labels.operations(1)), char(labels.operations(2))};
                bc1 = boxchart(ones(size(data_stim)), data_stim, 'BoxFaceColor', color_stim, 'BoxEdgeColor', 'k','LineWidth',1,'MarkerStyle','none','BoxFaceAlpha',1);hold on;
                bc2 = boxchart(2*ones(size(data_sham)), data_sham, 'BoxFaceColor', color_sham, 'BoxEdgeColor', 'k','LineWidth',1,'MarkerStyle','none','BoxFaceAlpha', 1);hold on;
                x_stim_scatter = 1 + (rand(size(data_stim))-0.5)*0.15;
                x_sham_scatter = 2 + (rand(size(data_sham))-0.5)*0.15;
                scatter(x_stim_scatter, data_stim, 10,'filled', ...
                    'MarkerFaceColor', [255,0,0]/255,'MarkerEdgeColor', [255,0,0]/255, 'LineWidth', 1, 'MarkerEdgeAlpha', 1);
                scatter(x_sham_scatter, data_sham, 10,'filled', ...
                    'MarkerFaceColor', [0,0,255]/255,'MarkerEdgeColor', [0,0,255]/255, 'LineWidth', 1, 'MarkerEdgeAlpha', 1);
                p_val = p_values(ch);title_str = sprintf('%s', labels.channels{ch});
                if p_val < 0.05
                    title(title_str, 'Color', 'r', 'FontWeight', 'bold');
                    sigstar({[1,2]}, p_val);
                end
                set(gca, 'XTick', [1, 2], 'XTickLabel', group_names, ...
                    'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
                xlim([0.5, 2.5]);title(labels.channels{ch}, 'Color', 'k');
                if mod(ch, grid_cols) == 1, ylabel('Power Change'); end
            else
                b = bar([1, 2], [m_stim,m_sham ], 'EdgeColor', 'flat', 'FaceColor', 'flat');hold on;
                set(b, 'FaceColor', 'none');
                color_val = [141,141,255]/255;
                b.CData(2,:) =[141,141,255]/255; 
                b.CData(1,:) = [255,128,128]/255; 
                b.EdgeColor = 'flat';
                b.LineWidth = 2;set(gca, 'TickDir', 'out', 'Box', 'off'); 
                x_sham_scatter = 2 + (rand(size(data_sham))-0.5)*0.3;
                x_stim_scatter = 1 + (rand(size(data_stim))-0.5)*0.3;
                scatter(x_sham_scatter, data_sham, 15, 'MarkerEdgeColor',[0,0,255]/255,'LineWidth',2);hold on;%'MarkerFaceAlpha', 0.3
                scatter(x_stim_scatter, data_stim, 15, 'MarkerEdgeColor',[255,0,0]/255,'LineWidth',2);hold on;%, 'filled', 'MarkerFaceAlpha', 0.3
                errorbar(2, m_sham, sem_sham,  'Color',[141,141,255]/255, 'LineWidth', 1.5, 'CapSize', 8);hold on;
                errorbar(1, m_stim, sem_stim, 'Color',[255,128,128]/255, 'LineWidth', 1.5, 'CapSize', 8);hold on;
                ch_title = labels.channels{ch};
                p_val = p_values(ch);
                title_str = sprintf('%s (p=%.3f)', ch_title, p_val);
                if p_val < 0.05
                    title(title_str, 'Color', 'r', 'FontWeight', 'bold');
                    sigstar({[1,2]}, p_val);
                else
                    title(title_str, 'Color', 'k');
                end
                xlim([0.5, 2.5]);
                set(gca, 'XTick', [1, 2], 'XTickLabel', labels.operations,'linewidth',2);
                box off;
            end
        end
        sgtitle(['Channel-wise Analysis for ' target_band_name], 'FontSize', 12, 'FontWeight', 'bold');
        exportgraphics(h, [labels.path,['Channel-wise Analysis for ' target_band_name],'.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
        saveas(h,[labels.path,['Channel-wise Analysis for ' target_band_name],'.png']);
        close(h);
    end
end

%% Correlation map
if length(behavior_stim1) ~= n_subs
    warning('Stim组 EEG 被试数 (%d) 与行为学数据 (%d) 不一致，请检查!', n_subs, length(behavior_stim));
end
for iii=1
    defaults.marker_color = [0, 0, 0];
    defaults.marker_size = 10;
    defaults.marker_face_alpha = 0.7;
    defaults.line_color = [0.8, 0.2, 0.2];
    defaults.line_width = 1;
    defaults.font_size = 5;
    defaults.grid_on = false;
    defaults.confidence_interval = true;
    defaults.show_equation = false;
    defaults.standardize = false;
    defaults.text_position = 'auto';
    corr_chan_name = {'Fp1','FC5','F3','T7','CP6'};
    for chan = 1:length(corr_chan_name)
        corr_band_idx = 3; % Alpha
        chan_idx = find(strcmpi(labels.channels, corr_chan_name{chan}));
        if isempty(chan_idx), error('未找到通道: %s', corr_chan_name{chan}); end

        eeg_stim_fz = Stim_Delta(:, corr_band_idx, chan_idx);
        eeg_sham_fz = Sham_Delta(:, corr_band_idx, chan_idx);

        X_stim = eeg_stim_fz;  Y_stim = behavior_stim;
        X_sham = eeg_sham_fz;  Y_sham = behavior_sham;
        X_all = [X_sham; X_stim]; Y_all = [Y_sham; Y_stim];

        defaults.xlabel = ['Δ Alpha (' corr_chan_name{chan} ')'];
        defaults.ylabel = 'Δ Arousal Score';

        h=figure('Position',[10,10,1500,400], 'Name', ['Corr: ' corr_chan_name{chan}]);
        subplot(1,3,1); plot_eeg_behavior_correlation(X_all, Y_all, defaults); title('All Subjects'); hold on;
        scatter(X_stim, Y_stim, 50, 'r', 'filled');
        scatter(X_sham, Y_sham, 50, 'k', 'filled');
        subplot(1,3,2); plot_eeg_behavior_correlation(X_stim, Y_stim, defaults); title([corr_chan_name{chan} ' Stim']);
        subplot(1,3,3); plot_eeg_behavior_correlation(X_sham, Y_sham, defaults); title([corr_chan_name{chan} ' Sham']);
        %saveas(h,[labels.path, ['Corr Alpha power' corr_chan_name{chan}],'.png']);close(h);
    end

    fprintf('\n正在进行组合通道脑-行为相关性分析...\n');
    analyze_bands = 1:11;power_bahav_stim={};power_bahav_sham={};
    num_total_bands = length(analyze_bands);
    num_rows = ceil(num_total_bands / 2); 
    for btype = 1:2
        figure('Color', 'w', 'Position', [50, 50, 1200, 150 * num_rows]);
        kk = 1; 
        for b_idx = 1:num_total_bands
            corr_band_idx = analyze_bands(b_idx);
            corr_chan_group = labels.channels(find(sig_mask(:,:,corr_band_idx)==1));
            if isempty(corr_chan_group), continue; end
            chan_idx_group = [];
            for ii=1:length(corr_chan_group)
                chan_idx_group(ii) = find(strcmpi(labels.channels, corr_chan_group{ii}));
            end
            eeg_stim_grp = squeeze(nansum(Stim_Delta(:, corr_band_idx, chan_idx_group), 3));
            eeg_sham_grp = squeeze(nansum(Sham_Delta(:, corr_band_idx, chan_idx_group), 3));
            X_stim = eeg_stim_grp; X_sham = eeg_sham_grp;
            switch btype
                case 1; behavior_stim=behavior_stim1; behavior_sham=behavior_sham1; labely='\Delta Arousal';
                case 2; behavior_stim=behavior_stim2; behavior_sham=behavior_sham2; labely='\Delta BDI';
            end
            Y_stim = behavior_stim; Y_sham = behavior_sham;
            X_all = [X_sham; X_stim]; Y_all = [Y_sham; Y_stim];
            defaults.title = labels.frequencies{corr_band_idx};defaults.ylabel = labely;defaults.xlabel = strcat(['Δ Power-',corr_chan_group{:}]);
            row_pos = ceil(b_idx / 2);
            col_offset = mod(b_idx-1, 2) * 3; 
            subplot(num_rows, 6, (row_pos-1)*6 + col_offset + 1);
            plot_eeg_behavior_correlation(X_all, Y_all, defaults); hold on;
            scatter(X_stim, Y_stim, defaults.marker_size, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
            scatter(X_sham, Y_sham, defaults.marker_size, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
            subplot(num_rows, 6, (row_pos-1)*6 + col_offset + 2);
            plot_eeg_behavior_correlation(X_stim, Y_stim, defaults);
            subplot(num_rows, 6, (row_pos-1)*6 + col_offset + 3);
            plot_eeg_behavior_correlation(X_sham, Y_sham, defaults);
            power_bahav_stim{btype,b_idx}=[X_stim, Y_stim];
            power_bahav_sham{btype,b_idx}=[X_sham, Y_sham];
        end
        set(gcf, 'Renderer', 'Painters');
    end
    save([labels.path,'power_bahav_corr.mat'],'power_bahav_stim','power_bahav_sham');
end
end

%% ==================== 辅助函数 ====================

function labels = generateLabels()
% (保持不变)
labels.path = 'C:\视频任务数据\论文数据整理\task2 data and code\Fig4_abc\Results\';
labels.operations = {'10Hz-Stim.','Sham'};
labels.operations_forload = {'10hz_psd_neg.mat','sham_psd_neg.mat'};
labels.frequencies = {['Delta'],['Theta'],['Alpha'],['LowBeta'],...
    ['HighBeta'],['Gamma1'],['Gamma2'],['Gamma3'],['Gamma4'],...
    ['Gamma5'],['Gamma6']};

labels.channels = {'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6',...
    'M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6',...
    'P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};
labels.region_mapping = {
    'Fp1', 'Frontal'; 'Fpz', 'Frontal'; 'Fp2', 'Frontal';
    'F7', 'Frontal'; 'F3', 'Frontal'; 'Fz', 'Frontal'; 'F4', 'Frontal'; 'F8', 'Frontal';
    'FC5', 'Fronto-Central'; 'FC1', 'Fronto-Central'; 'FC2', 'Fronto-Central'; 'FC6', 'Fronto-Central';
    'M1', 'Temporal'; 'T7', 'Temporal'; 'C3', 'Central'; 'Cz', 'Central'; 'C4', 'Central'; 'T8', 'Temporal'; 'M2', 'Temporal';
    'CP5', 'Centro-Parietal'; 'CP1', 'Centro-Parietal'; 'CP2', 'Centro-Parietal'; 'CP6', 'Centro-Parietal';
    'P7', 'Parietal'; 'P3', 'Parietal'; 'Pz', 'Parietal'; 'P4', 'Parietal'; 'P8', 'Parietal';
    'POz', 'Parieto-Occipital'; 'O1', 'Occipital'; 'Oz', 'Occipital'; 'O2', 'Occipital'};
unique_regions = unique(labels.region_mapping(:,2));
labels.regions = unique_regions';
labels.channel_regions = cell(size(labels.channels));
for i = 1:length(labels.channels)
    idx = find(strcmp(labels.region_mapping(:,1), labels.channels{i}));
    if ~isempty(idx)
        labels.channel_regions{i} = labels.region_mapping{idx,2};
    else
        labels.channel_regions{i} = 'Unknown';
    end
end
labels.colors.operations = [[255,0,0]; [0,0,255]]/255;
labels.colors.timepoints = [0.3 0.3 0.8; 0.8 0.3 0.3];
labels.colors.frequencies = [0.4  0.6  0.9];
labels.timepoints = {'Pre', 'Post'};
end