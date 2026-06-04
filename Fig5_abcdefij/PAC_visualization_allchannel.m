%% PAC-Behavior-Mechanism Integrated Analysis Script
clear; clc; close all;
% A=load('C:\视频任务数据\Network-index\PAC\通道间\所有通道-枕叶\hz10_pac_bottom_up.mat')
% B=load('C:\视频任务数据\Network-index\PAC\通道间\10hz_pac_z_m.mat','pac_z_score'); 
% legind=min(size(A.pac_z_score{1}{1},2),size(A.pac_z_score{2}{1},2));
% for indi=1:legind
%     for freqi=1:8
%         a=A.pac_z_score{1}{1}{indi}{freqi};size(a)
%         pac_z_score(1,indi,freqi,:,:,:,:)=a;
%         a=A.pac_z_score{2}{1}{indi}{freqi};size(a)
%         pac_z_score(2,indi,freqi,:,:,:,:)=a;
%     end
% end
% pac_z_score=cat(5,pac_z_score,B.pac_z_score(:,1:legind,:,7:9,:,:,:));
% chan_names = {'Fpz','Fp2','Fz','F4','F8','FC1','FC2','FC6','M1','C3','Cz','C4','T8','M2','CP5','CP1','CP2','P7','P3','Pz','P4','P8','POz','Fp1','F7','F3','FC5','T7','CP6'}; %
% labels.brain_regions = {'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6',...
%     'M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6',...
%     'P7','P3','Pz','P4','P8','POz'};
% chan_ind=[];
% for ii=1:length(labels.brain_regions)
%     chan_ind(ii)=find(strcmp(chan_names,labels.brain_regions{ii}));
% end
% pac_z_score=pac_z_score(:,:,:,:,chan_ind,:,:);
% save('C:\视频任务数据\Network-index\PAC\通道间\所有通道-枕叶\10hz_pac_z_m_bottom_up.mat','pac_z_score');

%% ==================== 1. 数据初始化与参数设置 ====================
isstim=1;
savepath='C:\视频任务数据\Network-index\PAC\通道间\所有通道-枕叶\topdown\';
if isstim==1;
    load('C:\视频任务数据\Network-index\PAC\通道间\所有通道-枕叶\10hz_pac_z_m_top_down.mat','pac_z_score'); 
elseif isstim==0;
    load('C:\视频任务数据\Network-index\PAC\通道间\所有通道-枕叶\sham_pac_z_m_top_down.mat','pac_z_score');
end
pac_z_score=squeeze(nanmean(pac_z_score,3));
n_sub = size(pac_z_score, 2);
alpha_axis = linspace(0.5,20,39);
gamma_axis = linspace(30,90,60);
load('C:\视频任务数据\videovalence.mat');
videovalence([[1,2,6],20+[1,5]],:) = [];
A = table2array(videovalence);
ii = 10; ind1 = find(A(:,2)==isstim);ind=ind1;
behavior_data = (A(ind,ii)-A(ind,ii-1))./(A(ind,ii-1));
load('C:\视频任务数据\BDI.mat');
BDI([[1,2,6,20],20+[1,5]],:) = []; 
A = table2array(BDI);
ii = 4; ind1 = find(A(:,2)==isstim);ind=ind1;
behavior_data2 = (A(ind,ii)-A(ind,ii-1))./(A(ind,ii-1)+0.001^1);
data_pre=squeeze(pac_z_score(1,:,:,:,:,:));data_post=squeeze(pac_z_score(2,:,:,:,:,:));size(data_post)
all_diff = (data_post - data_pre);
mean_diff_all = squeeze(nanmean(all_diff, 1));
chan_names = {'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6',...
    'M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6',...
    'P7','P3','Pz','P4','P8','POz'};
chan_names2={'O1','Oz','O2'} ;
n_chan = length(chan_names);
alpha_lvl = 0.05;
smooth_sigma = 1.2;
min_cluster_size=round(0.01*numel(alpha_axis)*numel(gamma_axis));
%% MI change
chanj=1:29;
chani=1:3;
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
sig_mask_alpha_gamma=[];
sig_mask_alpha_gamma_p={};
sig_mask_alpha_gamma_r={};
for chani=1;
    for chanj=1;
        a_sub = squeeze(nanmean(nanmean(data_pre(:,1:29,1:3, :, :), 3), 2));%[1,5,9,14,23]
        b_sub = squeeze(nanmean(nanmean(data_post(:,1:29,1:3, :, :), 3), 2));%[1,5,9,14,23]
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
            c_raw = [a(:); b(:)]; clims_raw = [min(c_raw), max(c_raw)];clims_raw=[-0.18 0.56];
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
            %% ==================== 5. 所有显著簇提取与行为相关性分析 ====================
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
                    case 2; behavior_datay=behavior_data2;labely='\Delta BDI';if isstim==1; sub_cluster_change(15,:)=[];end
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
                %close all;
            end
        end
    end
end
saveas(h,[savepath,'PAC change_all.png']);
exportgraphics(h, [savepath,'PAC change_all.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
%% the mean MI change across distinct cortical areas
% left frontal(Fpz, Fp1, AF3, AF7, Fz, F1, F3, F5, F7, FC1, FC3, FC5);
% right frontal (Fpz, Fp2, AF4, AF8, Fz, F2, F4, F6, F8, FC2, FC4, FC6);
% left central (Cz, C1, C3, C5, CPz, CP1, CP3, CP5) ;
% right central (Cz, C2, C4, C6, CPz, CP2, CP4, CP6);
% left temporal(FT7, T7, TP7) ;
% right temporal (FT8, T8, TP8);
% left occipital (Pz, P3, P5, P7, POz, PO3, PO7, Oz, O1) ;
% right occipital (Pz, P4, P6, P8, POz, PO4, PO8, Oz, O2);
regions.Prefrontal = 1:3;   % Fp1, Fpz, Fp2
regions.Frontal = 4:12;     % F7, F3, Fz, F4, F8, FC5, FC1, FC2, FC6
regions.Central = [14:18, 20:23]; % T7, C3, Cz, C4, T8, CP5, CP1, CP2, CP6
regions.Parietal = 24:28;   % P7, P3, Pz, P4, P8
regions.Occipital = 1:3;     % O1,O2,Oz !!!!!!!!!!!!!!!!!!!!是接受区域index
reg_names = fieldnames(regions);
num_regions = length(reg_names);
reg_sig_percentage = zeros(num_regions, num_regions);
reg_mean_change = zeros(num_regions, num_regions);
fprintf('正在统计各区域间的显著 PAC 变化率...\n');
for r_p = 1:num_regions 
    p_indices = regions.(reg_names{r_p});
    for r_a = 2:num_regions 
        try
            a_indices =regions.(reg_names{r_a});
            sig_count = 0; 
            total_pairs = length(p_indices) * length(a_indices);
            pair_changes = [];
            for i = p_indices
                for j = a_indices
                    curr_pre = squeeze(data_pre(:, i, j, :, :));
                    curr_post = squeeze(data_post(:, i, j, :, :));
                    p_map = zeros(size(curr_pre, 2), size(curr_pre, 3));
                    for fg = 1:size(curr_pre, 2)
                        for fa = 1:size(curr_pre, 3)
                            p_map(fg, fa) = signrank(curr_post(:, fg, fa), curr_pre(:, fg, fa));
                        end
                    end
                    [~, p_fdr]=mafdr(p_map(:));
                    p_fdr=reshape(p_fdr,size(p_map,1),size(p_map,2));
                    %binary_mask = p_fdr <= 0.05;
                    binary_mask = p_map < 0.05;
                    binary_mask = bwareaopen(binary_mask, min_cluster_size);
                    if any(binary_mask(:))
                        sig_count = sig_count + 1;
                    end
                    if true;%any(binary_mask(:))
                        diff_val = (curr_post - curr_pre);
                        pair_changes = [pair_changes, mean(diff_val(:))];
                    end
                end
            end
            reg_sig_percentage(r_p, r_a) = (sig_count / total_pairs) * 100;
            if ~isempty(pair_changes)
                reg_mean_change(r_p, r_a) = mean(pair_changes);
            end
        catch
            warning(strcat(reg_names{r_p},'-',reg_names{r_a},' not exist'))
            continue;
        end
    end
end
figure('Color', 'w', 'Name', 'Regional PAC Change Rate');
subplot(1,2,1);
imagesc(reg_sig_percentage);
colormap(hot); colorbar;
set(gca, 'XTick', 1:num_regions, 'XTickLabel', reg_names, ...
         'YTick', 1:num_regions, 'YTickLabel', reg_names);
title('Significant Connection Percentage (%)');
xlabel('Amplitude Regions'); ylabel('Phase Regions');
axis square;
subplot(1,2,2);
imagesc(reg_mean_change);
colormap(redblue); colorbar;
max_c = max(abs(reg_mean_change(:)));
caxis([-max_c, max_c]);
set(gca, 'XTick', 1:num_regions, 'XTickLabel', reg_names, ...
         'YTick', 1:num_regions, 'YTickLabel', reg_names);
title('Mean PAC Change Intensity (\Delta Z)');
xlabel('Amplitude Regions'); ylabel('Phase Regions');
axis square;
