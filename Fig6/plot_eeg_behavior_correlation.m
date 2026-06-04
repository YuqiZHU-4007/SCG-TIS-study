function [p_value,R_squared]=plot_eeg_behavior_correlation(eeg_index, behavior_improvement, options)
% 绘制脑电指数与行为改善相关性图（标准化版本）
% 输入：
%   eeg_index - 脑电指数数据（向量）
%   behavior_improvement - 行为改善数据（向量）
%   options - 可选参数结构体

% 参数检查
if nargin < 2
    error('至少需要两个输入参数：eeg_index和behavior_improvement');
end

% 数据验证
if ~isvector(eeg_index) || ~isvector(behavior_improvement)
    error('输入数据必须是向量');
end

if length(eeg_index) ~= length(behavior_improvement)
    error('两个输入向量的长度必须相同');
end

% 转换为列向量
eeg_index = eeg_index(:);
behavior_improvement = behavior_improvement(:);

% 移除NaN值
valid_idx = ~isnan(eeg_index) & ~isnan(behavior_improvement);
eeg_index = eeg_index(valid_idx);
behavior_improvement = behavior_improvement(valid_idx);

if length(eeg_index) < 3
    error('有效数据点少于3个，无法进行相关性分析');
end

% 数据标准化（z-score标准化）
eeg_index_z = normalize(eeg_index,'Range');
behavior_improvement_z = zscore(behavior_improvement);

% 设置默认参数
defaults.title = '';
defaults.xlabel = 'EEG（z-score）';
defaults.ylabel = 'Behavior（z-score）';
defaults.marker_color = [0.2, 0.4, 0.8];  % RGB颜色
defaults.marker_size = 60;
defaults.marker_face_alpha = 0.7;
defaults.line_color = [0.8, 0.2, 0.2];
defaults.line_width = 2;
defaults.font_size = 12;
defaults.grid_on = false;
defaults.confidence_interval = true;  % 是否显示置信区间
defaults.show_equation = false;  % 是否显示回归方程
defaults.text_position = 'auto';  % 'auto', 'top-right', 'top-left', 'bottom-right', 'bottom-left'
defaults.standardize = true;  % 是否标准化数据
defaults.show_pvalue=false;

% 合并用户提供的参数
if nargin < 3
    options = struct();
end
option_names = fieldnames(options);
for i = 1:length(option_names)
    defaults.(option_names{i}) = options.(option_names{i});
end
options = defaults;

% 选择使用标准化数据还是原始数据
if options.standardize
    x_data = eeg_index_z;
    y_data = behavior_improvement_z;
else
    x_data = eeg_index;
    y_data = behavior_improvement;
end

% 绘制散点图
h_scatter = scatter(x_data, y_data, options.marker_size, ...
    'MarkerFaceColor', options.marker_color, ...
    'MarkerEdgeColor', options.marker_color, ...
    'MarkerFaceAlpha', options.marker_face_alpha, ...
    'MarkerEdgeAlpha',0,...
    'LineWidth', 1.5, ...
    'DisplayName', '数据点');hold on;

% 计算线性回归
p = polyfit(x_data, y_data, 1);
y_fit = polyval(p, x_data);

% 计算R²和p值
[r, p_value] = corrcoef(x_data, y_data);
R_squared = r(1,2);%r(1,2)^2;
p_value = p_value(1,2);
% % 3. 执行 Spearman 相关分析 (非参数，对异常值更鲁棒)
%[R_squared, p_value] = corr(x_data, y_data, 'type', 'Spearman');

% 绘制回归线
x_range = linspace(min(x_data), max(x_data), 100);
y_fit_range = polyval(p, x_range);
h_line = plot(x_range, y_fit_range, 'Color', options.line_color, ...
    'LineWidth', options.line_width, ...
    'DisplayName', '回归线');

% 绘制曲线误差线（预测区间）
if options.confidence_interval
    % 计算预测区间
    alpha = 0.05;  % 95% 置信区间
    n = length(x_data);
    x_mean = mean(x_data);
    Sxx = sum((x_data - x_mean).^2);
    
    % 计算每个x点的标准误差
    residuals = y_data - y_fit;
    MSE = sum(residuals.^2) / (n - 2);
    
    % 预测区间的标准误差
    SE_pred = sqrt(MSE * (1 + 1/n + (x_range - x_mean).^2 / Sxx));
    
    % t分布的临界值
    t_critical = tinv(1 - alpha/2, n - 2);
    
    % 预测区间的上下限
    y_pred_upper = y_fit_range + t_critical * SE_pred;
    y_pred_lower = y_fit_range - t_critical * SE_pred;
    
    % 绘制曲线误差线（使用patch填充）
    h_fill = fill([x_range, fliplr(x_range)], ...
                  [y_pred_upper, fliplr(y_pred_lower)], ...
                  options.line_color, ...
                  'FaceAlpha', 0.15, ...
                  'EdgeColor', 'none', ...
                  'DisplayName', '95% 预测区间');
    
    % 绘制误差线的边界线（虚线）
    plot(x_range, y_pred_upper, '--', 'Color', options.line_color, ...
        'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot(x_range, y_pred_lower, '--', 'Color', options.line_color, ...
        'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% 计算文本标签的最佳位置（防重叠）
x_lim = xlim;
y_lim = ylim;
x_range_total = diff(x_lim);
y_range_total = diff(y_lim);

% 根据选项或自动选择文本位置
if strcmp(options.text_position, 'auto')
    % 自动选择数据点最稀疏的角落
    % 将图形分为4个象限
    x_mid = mean(x_lim);
    y_mid = mean(y_lim);
    
    % 计算每个象限的数据点密度
    density = zeros(1,4);  % [左上, 右上, 左下, 右下]
    
    % 左上象限
    density(1) = sum(x_data < x_mid & y_data > y_mid);
    % 右上象限
    density(2) = sum(x_data > x_mid & y_data > y_mid);
    % 左下象限
    density(3) = sum(x_data < x_mid & y_data < y_mid);
    % 右下象限
    density(4) = sum(x_data > x_mid & y_data < y_mid);
    
    [~, min_idx] = min(density);
    
    switch min_idx
        case 1  % 左上
            text_x = x_lim(1) + 0.05 * x_range_total;
            text_y = y_lim(2) - 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'left', 'VerticalAlignment', 'top'};
        case 2  % 右上
            text_x = x_lim(2) - 0.05 * x_range_total;
            text_y = y_lim(2) - 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'};
        case 3  % 左下
            text_x = x_lim(1) + 0.05 * x_range_total;
            text_y = y_lim(1) + 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom'};
        case 4  % 右下
            text_x = x_lim(2) - 0.05 * x_range_total;
            text_y = y_lim(1) + 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom'};
    end
else
    % 根据用户指定的位置
    switch options.text_position
        case 'top-right'
            text_x = x_lim(2) - 0.05 * x_range_total;
            text_y = y_lim(2) - 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'};
        case 'top-left'
            text_x = x_lim(1) + 0.05 * x_range_total;
            text_y = y_lim(2) - 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'left', 'VerticalAlignment', 'top'};
        case 'bottom-right'
            text_x = x_lim(2) - 0.05 * x_range_total;
            text_y = y_lim(1) + 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom'};
        case 'bottom-left'
            text_x = x_lim(1) + 0.05 * x_range_total;
            text_y = y_lim(1) + 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom'};
        otherwise
            text_x = x_lim(2) - 0.05 * x_range_total;
            text_y = y_lim(2) - 0.05 * y_range_total;
            text_align = {'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'};
    end
end

% 添加R²和p值文本（使用防重叠位置）
text_str = sprintf('R² = %.3f\np = %.4f', R_squared, p_value);labcolor='k';
if p_value < 0.05
    text_str = sprintf('R² = %.3f\np = %.4f', R_squared, p_value);%text_str = sprintf('R² = %.3f\np < 0.05', R_squared);
    labcolor='r';
end
if options.show_pvalue
h_text = text(text_x, text_y, text_str, ...
    'FontSize', options.font_size + 2, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', [1, 1, 1, 0.85], ...
    'EdgeColor', [0.3, 0.3, 0.3], ...
    'LineWidth', 1, ...
    'Margin', 8, ...
    text_align{:});
end
% 添加回归方程（放在统计文本下方）
if options.show_equation
    eq_text = sprintf('y = %.3fx + %.3f', p(1), p(2));
    
    % 根据文本对齐方式确定方程位置
    switch text_align{2}
        case 'right'
            eq_align = {'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'};
            eq_x = text_x;
            eq_y = text_y - 0.08 * y_range_total;
        case 'left'
            eq_align = {'HorizontalAlignment', 'left', 'VerticalAlignment', 'top'};
            eq_x = text_x;
            eq_y = text_y - 0.08 * y_range_total;
    end
    
    h_eq = text(eq_x, eq_y, eq_text, ...
        'FontSize', options.font_size, ...
        'BackgroundColor', [1, 1, 1, 0.85], ...
        'EdgeColor', [0.3, 0.3, 0.3], ...
        'LineWidth', 0.5, ...
        'Margin', 4, ...
        eq_align{:});
end

% 设置图形属性
xlabel(options.xlabel, 'FontSize', options.font_size + 1, 'FontWeight', 'bold');
ylabel(options.ylabel, 'FontSize', options.font_size + 1, 'FontWeight', 'bold');
title([options.title,'-',text_str], 'FontSize', options.font_size + 3, 'FontWeight', 'bold','Color',labcolor);

if options.grid_on
    grid on;
    grid minor;
    set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
end

% 添加图例
% if options.confidence_interval
%     legend([h_scatter, h_line, h_fill], 'Location', 'best', ...
%         'FontSize', options.font_size - 1);
% else
%     legend([h_scatter, h_line], 'Location', 'best', ...
%         'FontSize', options.font_size - 1);
% end

box off;
hold off;

% 设置坐标轴字体和美化
set(gca, 'FontSize', options.font_size, ...
    'TickDir', 'out', ...
    'Box', 'on', ...
    'LineWidth', 1);

% 设置图形背景色为白色
set(gcf, 'Color', 'w');

% 调整坐标轴范围，为文本留出空间
padding_x = 0.08 * diff(x_lim);
padding_y = 0.08 * diff(y_lim);
xlim([x_lim(1)-padding_x, x_lim(2)+padding_x]);
ylim([y_lim(1)-padding_y, y_lim(2)+padding_y]);

% % 输出统计结果到命令行
% fprintf('\n========== 相关性分析结果 ==========\n');
% fprintf('数据点数: %d\n', length(x_data));
% if options.standardize
%     fprintf('数据已标准化\n');
% end
% fprintf('相关系数 r: %.4f\n', r(1,2));
% fprintf('决定系数 R²: %.4f\n', R_squared);
% fprintf('p值: %.6f\n', p_value);
% fprintf('回归方程: y = %.4f*x + %.4f\n', p(1), p(2));
% fprintf('====================================\n\n');

end
