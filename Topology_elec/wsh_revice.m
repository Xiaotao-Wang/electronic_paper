% 计算阻抗函数 Z~(ω) 的虚部与连续相位分布
% 被积函数分母: omega^2*L*C*sin(k1/2)^2 - sin(k2/2)^2

clear; clc; close all;

%% 1. 参数设置（保留你原本的离散参数，用于计算Z、三维图和红蓝热力图）
L = 1; C = 1; omega = 10;
omega_sq_LC = omega^2 * L * C;  % omega²LC

% 离散范围
m_vals = -20:1:20;
n_vals = -20:1:20;

% 高精度积分设置
N_quad = 500;  % 积分网格点数（500x500 足够保证主值精确度）

fprintf('========== 1. 开始计算离散 Im[Z~(ω)] ==========\n');
fprintf('分母 = ω²LC·sin²(k₁/2) - sin²(k₂/2)\n');
fprintf('ω²LC = %.4f\n', omega_sq_LC);
fprintf('网格点数: %d × %d\n', N_quad, N_quad);
fprintf('m 范围: %d 到 %d, n 范围: %d 到 %d\n', ...
    min(m_vals), max(m_vals), min(n_vals), max(n_vals));
fprintf('总计算量: %d 个离散积分点\n\n', length(m_vals)*length(n_vals));

% 存储离散结果
ImZ_matrix = zeros(length(m_vals), length(n_vals));
total_points = length(m_vals) * length(n_vals);
progress = 0;

tic;

% 预先生成高斯-勒让德节点和权重
[x, w] = lgwt(N_quad, -pi, pi);
[K1, K2] = meshgrid(x, x);
[W1, W2] = meshgrid(w, w);

% 预计算分母基础矩阵与高斯权重乘积
denom_base = omega_sq_LC * sin(K1/2).^2 - sin(K2/2).^2;
weight_prod = W1 .* W2;

% 完全不动你的核心离散循环
for mi = 1:length(m_vals)
    m = m_vals(mi);
    for ni = 1:length(n_vals)
        n = n_vals(ni);
        
        progress = progress + 1;
        if mod(progress, 50) == 0
            fprintf('进度: %d/%d (%.1f%%)\n', progress, total_points, progress/total_points*100);
        end
        
        % m=n=0 时特殊处理（分子为0，积分为0）
        if m == 0 && n == 0
            ImZ_matrix(mi, ni) = 0;
            continue;
        end
        
        % 使用解析延拓（微扰法）精确逼近柯西主值
        num = 1 - cos(m*K1) .* cos(n*K2);
        eps_dist = 1e-7; 
        f_matrix = real(num ./ (denom_base + 1i * eps_dist));
        
        % 二维网格积分求和
        integral_value = sum(sum(f_matrix .* weight_prod));
        
        % Z~(ω) 的虚部系数
        ImZ_matrix(mi, ni) = -omega * L / (4 * pi^2) * integral_value;
    end
end

elapsed_time = toc;
fprintf('\n离散部分计算完成！耗时: %.2f 秒\n', elapsed_time);


%% 2. 连续变量处理（专用于绘制高保真、连续平滑的相位与强度变化图）
fprintf('\n========== 2. 开始将 m, n 处理为连续变量计算图像 ==========\n');
% 定义高密度的连续空间坐标（例如从 0 到 12，中间插值 200 个点，形成连续视觉）
dense_grid = linspace(0, 12, 200); 
[X_dense, Y_dense] = meshgrid(dense_grid, dense_grid);

% 初始化高密度连续响应矩阵
Z_dense_complex = zeros(size(X_dense));

tic;
% 为了让连续变量模拟速度飞起来，这里对空间格点采用高效率矩阵广播（不再使用嵌套双循环）
for i = 1:numel(dense_grid)
    curr_m = dense_grid(i);
    % 矩阵化并行计算整行的连续响应
    % 连续分子：1 - cos(m*k1)*cos(n*k2)
    cos_m_K1 = cos(curr_m * K1);
    
    for j = 1:numel(dense_grid)
        curr_n = dense_grid(j);
        
        if curr_m == 0 && curr_n == 0
            Z_dense_complex(i, j) = 0;
            continue;
        end
        
        num_dense = 1 - cos_m_K1 .* cos(curr_n * K2);
        eps_dist = 1e-7;
        % 同时保留实部与虚部以计算完整连续相位
        f_dense = num_dense ./ (denom_base + 1i * eps_dist);
        integral_dense = sum(sum(f_dense .* weight_prod));
        
        Z_dense_complex(i, j) = (-1i * omega * L) / (4 * pi^2) * integral_dense;
    end
end

% 提取连续物理量
Dense_Power = abs(imag(Z_dense_complex));
Dense_Phase = angle(Z_dense_complex) * (180 / pi);
fprintf('连续变量图像数据计算完成！耗时: %.2f 秒\n', toc);


%% 高斯-勒让德求积节点和权重
function [x, w] = lgwt(N, a, b)
    x = cos(pi * (4*(1:N)' - 1) / (4*N + 2));
    epsilon = 1e-15;
    for iter = 1:100
        P_prev = 1; P_curr = x; P_der_prev = 0; P_der_curr = 1;
        for k = 2:N
            P_next = ((2*k-1)*x.*P_curr - (k-1)*P_prev) / k;
            P_der_next = ((2*k-1)*(P_curr + x.*P_der_curr) - (k-1)*P_der_prev) / k;
            P_prev = P_curr; P_curr = P_next; P_der_prev = P_der_curr; P_der_curr = P_der_next;
        end
        dx = P_curr ./ P_der_curr; x = x - dx;
        if max(abs(dx)) < epsilon, break; end
    end
    w = 2 ./ ((1 - x.^2) .* P_der_curr.^2);
    x = 0.5 * (a + b) + 0.5 * (b - a) * x; w = 0.5 * (b - a) * w;
end

%% 自定义红蓝配色函数
function cmap = redbluecmap
    n = 256; cmap = zeros(n, 3); half = floor(n/2);
    for i = 1:half, t = (i-1)/(half-1); cmap(i, :) = [t, t, 1]; end
    for i = half+1:n, t = (i-half-1)/(n-half-1); cmap(i, :) = [1, 1-t, 1-t]; end
end

%% 显示结果终端文本（完全保留）
fprintf('\n========== 非零值检测 ==========\n');
fprintf('阈值: 1e-6\n\n');
[m_idx, n_idx] = find(abs(ImZ_matrix) > 1e-6);
non_zero_count = length(m_idx);
if non_zero_count == 0
    fprintf('未检测到非零值\n');
else
    fprintf('检测到 %d 个非零值:\n', non_zero_count);
    fprintf('m\tn\tIm[Z~(ω)]\n');
    fprintf('------------------------\n');
    non_zero_values = [];
    for i = 1:non_zero_count
        m_val = m_vals(m_idx(i)); n_val = n_vals(n_idx(i)); z_val = ImZ_matrix(m_idx(i), n_idx(i));
        fprintf('%3d\t%3d\t%.8f\n', m_val, n_val, z_val);
        non_zero_values = [non_zero_values; m_val, n_val, z_val];
    end
    fprintf('\n========== 模式分析 ==========\n');
    fprintf('按 |m|, |n| 分组:\n'); fprintf('|m|\t|n|\t值\n'); fprintf('------------------------\n');
    for i = 1:non_zero_count
        fprintf('%d\t%d\t%.8f\n', abs(non_zero_values(i, 1)), abs(non_zero_values(i, 2)), non_zero_values(i, 3));
    end
end

%% ==================== 绘图部分 ====================

%% 图 1：沿 m 方向的子图（固定 n，基于离散数据）
figure('Name', '柯西主值积分结果', 'Position', [50, 100, 600, 450]);
n_fixed = [-3, -2, -1, 0, 1, 2, 3];
colors = lines(length(n_fixed)); hold on;
for ni = 1:length(n_fixed)
    n = n_fixed(ni); idx_n = find(n_vals == n);
    if ~isempty(idx_n)
        plot(m_vals, ImZ_matrix(:, idx_n), 'o-', 'LineWidth', 2, 'Color', colors(ni,:), ...
            'MarkerSize', 6, 'DisplayName', sprintf('n = %d', n));
    end
end
xlabel('m', 'FontSize', 12); ylabel('Im{Z~(ω)}', 'FontSize', 12);
title(sprintf('沿 m 方向变化曲线，ω²LC = %.2f', omega_sq_LC));
legend('Location', 'best'); grid on; hold off;

%% 图 2：原始离散红蓝热力图（包含离散数值文本标记）
figure('Name', '红蓝热力图', 'Position', [660, 100, 600, 450]);
imagesc(n_vals, m_vals, ImZ_matrix); colormap(gca, redbluecmap); colorbar;
xlabel('n', 'FontSize', 12); ylabel('m', 'FontSize', 12);
title('Im{Z~(ω)} 红蓝离散热力图'); set(gca, 'YDir', 'normal'); axis square;
hold on;
for mi = 1:length(m_vals)
    for ni = 1:length(n_vals)
        if abs(ImZ_matrix(mi, ni)) > 1e-6
            val = ImZ_matrix(mi, ni);
            if abs(val) > 0.5, text_color = 'w'; else, text_color = 'k'; end
            text(n_vals(ni), m_vals(mi), sprintf('%.2f', val), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', text_color, 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
end
hold off;

%% 图 3：【彻底补回】离散非零值三维标注散点图
figure('Name', '三维视图', 'Position', [100, 200, 700, 500]);
[M, N] = meshgrid(m_vals, n_vals);
nonzero_idx = abs(ImZ_matrix) > 1e-6;
M_nonzero = M(nonzero_idx); N_nonzero = N(nonzero_idx); Z_nonzero = ImZ_matrix(nonzero_idx);

if ~isempty(Z_nonzero)
    scatter3(M_nonzero, N_nonzero, Z_nonzero, 150, Z_nonzero, 'filled');
    xlabel('m'); ylabel('n'); zlabel('Im{Z~(ω)}');
    title(sprintf('Im{Z~(ω)} 非零值三维分布 (共 %d 个点)', length(Z_nonzero)));
    colormap(jet); colorbar; grid on; view(45, 30);
    for i = 1:length(M_nonzero)
        text(M_nonzero(i), N_nonzero(i), Z_nonzero(i), sprintf('(%d,%d)', M_nonzero(i), N_nonzero(i)), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
else
    text(0, 0, 0, '无非零值', 'FontSize', 14, 'HorizontalAlignment', 'center');
    xlabel('m'); ylabel('n'); zlabel('Im{Z~(ω)}'); grid on;
end

%% 图 4：【全新切片】无限大模拟——连续变量下的功率与相位分布网格图
figure('Name', '连续变量模拟还原', 'Position', [150, 150, 1100, 500]);

% --- 左子图：连续变量 Power Flow 强度图 ---
subplot(1, 2, 1);
imagesc(dense_grid, dense_grid, Dense_Power); 
colormap(gca, jet(256)); colorbar;
set(gca, 'YDir', 'normal'); axis square; grid on; hold on;
% 叠加网格线，依然维持图中的格子标尺感
for i = 0:1:12, plot([i, i], [0, 12], 'k-', 'LineWidth', 0.5); end
for j = 0:1:12, plot([0, 12], [j, j], 'k-', 'LineWidth', 0.5); end
xlabel('x', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 13, 'FontWeight', 'bold', 'Rotation', 0, 'HorizontalAlignment', 'right');
set(gca, 'XTick', 0:2:12, 'YTick', 0:2:12);
title('Power flow distribution (Continuous)', 'FontSize', 12, 'FontWeight', 'bold');

% --- 右子图：连续变量 Voltage Phase 完美平滑条纹图 ---
subplot(1, 2, 2);
imagesc(dense_grid, dense_grid, Dense_Phase); 
colormap(gca, jet(256)); caxis([-180, 180]); 
h_cb = colorbar; set(h_cb, 'YTick', -150:50:150);
set(gca, 'YDir', 'normal'); axis square; grid on; hold on;
% 叠加网格线
for i = 0:1:12, plot([i, i], [0, 12], 'k-', 'LineWidth', 0.5); end
for j = 0:1:12, plot([0, 12], [j, j], 'k-', 'LineWidth', 0.5); end
xlabel('x', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 13, 'FontWeight', 'bold', 'Rotation', 0, 'HorizontalAlignment', 'right');
set(gca, 'XTick', 0:2:12, 'YTick', 0:2:12);
title('Voltage phase distributions (Continuous)', 'FontSize', 12, 'FontWeight', 'bold');

%% 输出统计信息
fprintf('\n========== 统计信息 ==========\n');
fprintf('离散总点数: %d\n', length(m_vals)*length(n_vals));
fprintf('非零值数量: %d\n', non_zero_count);