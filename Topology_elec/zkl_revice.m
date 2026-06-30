% 计算固定 (m,n) 时，积分值随 omega^2*L*C 的变化
% 使用高精度柯西主值积分

clear; clc; close all;

%% 1. 参数设置
L = 1; C = 1;
m_fixed = 1;  % 固定 m
n_fixed = 1;  % 固定 n

% omega^2*L*C 的变化范围
omega_sq_LC_vals = linspace(0.95, 1.05, 1000);

% 高精度积分设置
N_quad = 1000;  % 积分网格点数

fprintf('========== 高精度柯西主值积分计算 ==========\n');
fprintf('固定 (m=%d, n=%d) 积分值随 ω²LC 的变化\n', m_fixed, n_fixed);
fprintf('分母 = ω²LC·sin²(k₁/2) - sin²(k₂/2)\n');
fprintf('积分网格: %d × %d (高斯-勒让德求积)\n', N_quad, N_quad);

% 生成高斯-勒让德节点和权重
[x, w] = lgwt(N_quad, -pi, pi);
[K1, K2] = meshgrid(x, x);
[W1, W2] = meshgrid(w, w);
weight_prod = W1 .* W2;

% 定义分子（与 ω²LC 无关）
num = 1 - cos(m_fixed*K1) .* cos(n_fixed*K2);

% 存储结果
integral_values = zeros(size(omega_sq_LC_vals));
ImZ_values = zeros(size(omega_sq_LC_vals));

tic;

%% 对每个 omega^2*L*C 值计算柯西主值积分
for idx = 1:length(omega_sq_LC_vals)
    omega_sq_LC = omega_sq_LC_vals(idx);
    
    % 计算分母
    denom = omega_sq_LC * sin(K1/2).^2 - sin(K2/2).^2;
    
    % 使用虚部微扰计算主值
    eps_dist = 1e-8;
    f_matrix = real(num ./ (denom + 1i * eps_dist));
    
    % 检测并处理奇异点
    denom_abs = abs(denom);
    singular_threshold = 1e-3;
    singular_mask = denom_abs < singular_threshold;
    
    if any(singular_mask(:))
        [sing_i, sing_j] = find(singular_mask);
        
        for s = 1:min(length(sing_i), 20)
            i0 = sing_i(s);
            j0 = sing_j(s);
            
            radius = 8;
            i_range = max(1, i0-radius):min(N_quad, i0+radius);
            j_range = max(1, j0-radius):min(N_quad, j0+radius);
            
            if length(i_range) > 2 && length(j_range) > 2
                K1_local = K1(i_range, j_range);
                K2_local = K2(i_range, j_range);
                denom_local = omega_sq_LC * sin(K1_local/2).^2 - sin(K2_local/2).^2;
                num_local = 1 - cos(m_fixed*K1_local) .* cos(n_fixed*K2_local);
                
                valid_mask = abs(denom_local) > 1e-4;
                if sum(valid_mask(:)) > 10
                    f_local = real(num_local ./ (denom_local + 1i * 1e-6));
                    f_local(~valid_mask) = mean(f_local(valid_mask));
                    f_matrix(i_range, j_range) = f_local;
                end
            end
        end
    end
    
    % 二维网格积分求和
    integral_value = sum(sum(f_matrix .* weight_prod));
    
    % 检查异常值
    if abs(integral_value) > 100
        eps_dist = 1e-4;
        f_matrix = real(num ./ (denom + 1i * eps_dist));
        integral_value = sum(sum(f_matrix .* weight_prod));
    end
    
    integral_values(idx) = integral_value;
    
    % Z~(ω) 的虚部系数
    ImZ_values(idx) = -omega_sq_LC / (4 * pi^2) * integral_value;
    
    % 显示进度
    if mod(idx, 10) == 0
        fprintf('进度: %d/%d (%.1f%%), 当前值: %.4f\n', ...
            idx, length(omega_sq_LC_vals), idx/length(omega_sq_LC_vals)*100, ImZ_values(idx));
    end
end

elapsed_time = toc;
fprintf('\n计算完成！耗时: %.2f 秒\n', elapsed_time);

%% ==================== 绘图部分 ====================

%% 图 1：沿 omega_sq_LC 方向的子图
figure('Name', '柯西主值积分结果', 'Position', [50, 100, 1000, 450]);

% 左子图：积分值随 ω²LC 变化
subplot(1, 2, 1);
plot(omega_sq_LC_vals, integral_values, 'b-', 'LineWidth', 2);
hold on;
plot(omega_sq_LC_vals, zeros(size(omega_sq_LC_vals)), 'k--', 'LineWidth', 1);
xlabel('ω²LC', 'FontSize', 12);
ylabel('积分值 (主值)', 'FontSize', 12);
title(sprintf('柯西主值积分随 ω²LC 变化 (m=%d, n=%d)', m_fixed, n_fixed));
grid on;
legend('积分值', '零点', 'Location', 'best');
hold off;

% 右子图：Im[Z~(ω)] 随 ω²LC 变化
subplot(1, 2, 2);
plot(omega_sq_LC_vals, ImZ_values, 'r-', 'LineWidth', 2);
xlabel('ω²LC', 'FontSize', 12);
ylabel('Im{Z~(ω)}', 'FontSize', 12);
title(sprintf('Im{Z~(ω)} 随 ω²LC 变化 (m=%d, n=%d)', m_fixed, n_fixed));
grid on;
% 显示数值范围
ylim_vals = ylim;
text(0.05*max(omega_sq_LC_vals), 0.8*ylim_vals(2), ...
    sprintf('最大值: %.2f\n最小值: %.2f', max(ImZ_values), min(ImZ_values)), ...
    'FontSize', 10, 'BackgroundColor', 'w');

%% 图 2：红蓝热力图
figure('Name', '红蓝热力图', 'Position', [660, 100, 600, 450]);

% 构建二维矩阵
m_range = -3:1:3;
ImZ_matrix = zeros(length(m_range), length(omega_sq_LC_vals));

fprintf('\n计算热力图数据...\n');
for mi = 1:length(m_range)
    m = m_range(mi);
    num_temp = 1 - cos(m*K1) .* cos(n_fixed*K2);
    
    for idx = 1:length(omega_sq_LC_vals)
        omega_sq_LC = omega_sq_LC_vals(idx);
        denom = omega_sq_LC * sin(K1/2).^2 - sin(K2/2).^2;
        
        eps_dist = 1e-6;
        f_matrix = real(num_temp ./ (denom + 1i * eps_dist));
        
        denom_abs = abs(denom);
        singular_mask = denom_abs < 1e-3;
        if any(singular_mask(:))
            [sing_i, sing_j] = find(singular_mask);
            for s = 1:min(length(sing_i), 20)
                i0 = sing_i(s);
                j0 = sing_j(s);
                radius = 5;
                i_range = max(1, i0-radius):min(N_quad, i0+radius);
                j_range = max(1, j0-radius):min(N_quad, j0+radius);
                if length(i_range) > 2 && length(j_range) > 2
                    denom_local = omega_sq_LC * sin(K1(i_range, j_range)/2).^2 - sin(K2(i_range, j_range)/2).^2;
                    num_local = 1 - cos(m*K1(i_range, j_range)) .* cos(n_fixed*K2(i_range, j_range));
                    f_local = real(num_local ./ (denom_local + 1i * 1e-6));
                    valid_mask = abs(denom_local) > 1e-4;
                    if sum(valid_mask(:)) > 5
                        f_local(~valid_mask) = mean(f_local(valid_mask));
                        f_matrix(i_range, j_range) = f_local;
                    end
                end
            end
        end
        
        integral_value = sum(sum(f_matrix .* weight_prod));
        ImZ_matrix(mi, idx) = -omega_sq_LC / (4 * pi^2) * integral_value;
    end
    if mod(mi, 2) == 0
        fprintf('  m = %d 完成\n', m);
    end
end

% 限制颜色范围避免异常值影响显示
clim_max = max(abs(ImZ_matrix(:)));
if clim_max > 10
    clim_max = 10;
end

imagesc(omega_sq_LC_vals, m_range, ImZ_matrix);
colormap(gca, redbluecmap);
caxis([-clim_max, clim_max]);
colorbar;
xlabel('ω²LC', 'FontSize', 12);
ylabel('m', 'FontSize', 12);
title(sprintf('Im{Z~(ω)} 红蓝热力图 (固定 n=%d)', n_fixed));
set(gca, 'YDir', 'normal');
axis square;

%% 图 3：三维散点图
figure('Name', '三维视图', 'Position', [100, 200, 700, 500]);

[M_grid, W_grid] = meshgrid(m_range, omega_sq_LC_vals);
ImZ_flat = ImZ_matrix';
nonzero_idx = abs(ImZ_flat) > 1e-3;
M_nonzero = M_grid(nonzero_idx);
W_nonzero = W_grid(nonzero_idx);
Z_nonzero = ImZ_flat(nonzero_idx);

if ~isempty(Z_nonzero)
    scatter3(M_nonzero, W_nonzero, Z_nonzero, 100, Z_nonzero, 'filled');
    xlabel('m');
    ylabel('ω²LC');
    zlabel('Im{Z~(ω)}');
    title(sprintf('Im{Z~(ω)} 三维分布 (固定 n=%d)', n_fixed));
    colormap(jet);
    colorbar;
    grid on;
    view(45, 30);
else
    text(0, 0, 0, '无非零值', 'FontSize', 14, 'HorizontalAlignment', 'center');
    xlabel('m');
    ylabel('ω²LC');
    zlabel('Im{Z~(ω)}');
    grid on;
end

%% 图 4：连续变量模拟
figure('Name', '连续变量模拟', 'Position', [150, 150, 1100, 500]);

% 左子图：强度图
subplot(1, 2, 1);
imagesc(omega_sq_LC_vals, m_range, ImZ_matrix);
colormap(gca, jet(256));
caxis([-clim_max, clim_max]);
colorbar;
set(gca, 'YDir', 'normal');
axis square;
grid on;
xlabel('ω²LC', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('m', 'FontSize', 13, 'FontWeight', 'bold');
title('Im{Z~(ω)} 强度分布', 'FontSize', 12, 'FontWeight', 'bold');

% 右子图：相位分布
subplot(1, 2, 2);
Phase_matrix = angle(ImZ_matrix + 1i*eps) * (180/pi);
imagesc(omega_sq_LC_vals, m_range, Phase_matrix);
colormap(gca, jet(256));
caxis([-180, 180]);
h_cb = colorbar;
set(h_cb, 'YTick', -150:50:150);
set(gca, 'YDir', 'normal');
axis square;
grid on;
xlabel('ω²LC', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('m', 'FontSize', 13, 'FontWeight', 'bold');
title('相位分布', 'FontSize', 12, 'FontWeight', 'bold');

%% 输出统计信息
fprintf('\n========== 统计信息 ==========\n');
fprintf('ω²LC 总点数: %d\n', length(omega_sq_LC_vals));
fprintf('积分网格: %d × %d\n', N_quad, N_quad);
fprintf('Im[Z] 范围: %.4f 到 %.4f\n', min(ImZ_values), max(ImZ_values));
fprintf('Im[Z] 均值: %.4f\n', mean(ImZ_values));
fprintf('Im[Z] 标准差: %.4f\n', std(ImZ_values));

%% 自定义红蓝配色函数
function cmap = redbluecmap
    n = 256;
    cmap = zeros(n, 3);
    half = floor(n/2);
    for i = 1:half
        t = (i-1)/(half-1);
        cmap(i, :) = [t, t, 1];
    end
    for i = half+1:n
        t = (i-half-1)/(n-half-1);
        cmap(i, :) = [1, 1-t, 1-t];
    end
end

%% 高斯-勒让德求积节点和权重
function [x, w] = lgwt(N, a, b)
    x = cos(pi * (4*(1:N)' - 1) / (4*N + 2));
    epsilon = 1e-15;
    for iter = 1:100
        P_prev = 1;
        P_curr = x;
        P_der_prev = 0;
        P_der_curr = 1;
        for k = 2:N
            P_next = ((2*k-1)*x.*P_curr - (k-1)*P_prev) / k;
            P_der_next = ((2*k-1)*(P_curr + x.*P_der_curr) - (k-1)*P_der_prev) / k;
            P_prev = P_curr;
            P_curr = P_next;
            P_der_prev = P_der_curr;
            P_der_curr = P_der_next;
        end
        dx = P_curr ./ P_der_curr;
        x = x - dx;
        if max(abs(dx)) < epsilon
            break;
        end
    end
    w = 2 ./ ((1 - x.^2) .* P_der_curr.^2);
    x = 0.5 * (a + b) + 0.5 * (b - a) * x;
    w = 0.5 * (b - a) * w;
end