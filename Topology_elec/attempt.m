%% 计算固定 (m,n) 时，积分值随 omega^2*L*C 的变化
% 正确区分 ω²LC=1 和 ω²LC≠1 两种情况

clear; clc; close all;

%% 1. 参数设置
L = 1; C = 1;
m_fixed = 0;  % 固定 m
n_fixed = 5;  % 固定 n

% omega^2*L*C 的变化范围（避免恰好等于1，单独处理）
omega_sq_LC_vals = [linspace(0.999, 0.9999, 500), 1.0, linspace(1.0001, 1.001, 500)];

% 积分设置
N_quad = 400;  % 每个方向的积分点数

fprintf('========== 柯西主值积分计算（修正版）==========\n');
fprintf('固定 (m=%d, n=%d)\n', m_fixed, n_fixed);
fprintf('特别处理 ω²LC = 1 的情况\n\n');

% 生成均匀积分网格
k1_vec = linspace(-pi, pi, N_quad);
k2_vec = linspace(-pi, pi, N_quad);
dk1 = k1_vec(2) - k1_vec(1);
dk2 = k2_vec(2) - k2_vec(1);
[K1, K2] = meshgrid(k1_vec, k2_vec);
dA = dk1 * dk2;

% 定义分子
num = 1 - cos(m_fixed*K1) .* cos(n_fixed*K2);

% 存储结果
integral_values = zeros(size(omega_sq_LC_vals));
ImZ_values = zeros(size(omega_sq_LC_vals));
method_used = cell(size(omega_sq_LC_vals));

fprintf('开始计算...\n');
tic;

%% 对每个 omega^2*L*C 值计算柯西主值积分
for idx = 1:length(omega_sq_LC_vals)
    omega_sq_LC = omega_sq_LC_vals(idx);
    
    % ========== 情况1：ω²LC = 1（特殊处理）==========
    if abs(omega_sq_LC - 1.0) < 1e-15
        fprintf('========== 处理 ω²LC = 1.000000 ==========\n');
        method_used{idx} = 'ω²LC=1特殊处理';
        
        % 此时分母为：sin²(k₁/2) - sin²(k₂/2)
        % 奇点在 k₁ = ±k₂ + 2πn 处
        denom = sin(K1/2).^2 - sin(K2/2).^2;
        
        % 分析奇点结构
        fprintf('分母范围: [%.2e, %.2e]\n', min(denom(:)), max(denom(:)));
        
        % 方法1：分析积分 - 利用对称性和解析性质
        % 由于被积函数在四个象限具有对称性，可以将积分区域缩小到第一象限
        
        % 先尝试去除奇点线的数值积分
        % 奇点线: k₁ = k₂ 和 k₁ = -k₂
        % 设定一个去除奇点线的宽度 δ
        
        delta_vals = logspace(-6, -1, 100);  % 不同的去除宽度
        pv_estimates = zeros(size(delta_vals));
        num_excluded = zeros(size(delta_vals));
        
        for di = 1:length(delta_vals)
            delta = delta_vals(di);
            
            % 识别奇点线附近的点
            % 奇点条件: |k₁ - k₂| < δ 或 |k₁ + k₂| < δ
            % 同时考虑周期性: |k₁ - k₂ - 2π| < δ 等
            
            % 主对角线: k₁ ≈ k₂
            near_main_diag = abs(K1 - K2) < delta | abs(K1 - K2 - 2*pi) < delta | abs(K1 - K2 + 2*pi) < delta;
            % 反对角线: k₁ ≈ -k₂
            near_anti_diag = abs(K1 + K2) < delta | abs(K1 + K2 - 2*pi) < delta | abs(K1 + K2 + 2*pi) < delta;
            
            exclude_mask = near_main_diag | near_anti_diag;
            num_excluded(di) = sum(exclude_mask(:));
            
            % 如果排除的点太多（超过40%），这个delta太大，不适用
            if num_excluded(di) > 0.4 * N_quad^2
                pv_estimates(di) = NaN;
                continue;
            end
            
            % 计算排除奇点区域后的积分
            valid_mask = ~exclude_mask;
            f_valid = num(valid_mask) ./ denom(valid_mask);
            pv_estimates(di) = sum(f_valid(:)) * dA;
        end
        
        % 分析收敛性
        valid_idx = ~isnan(pv_estimates);
        if sum(valid_idx) > 5
            % 取delta较小但排除点不太多的区间
            % 选择num_excluded在合理范围内（1%-20%）的估计值
            reasonable_idx = valid_idx & (num_excluded > 0.01*N_quad^2) & (num_excluded < 0.2*N_quad^2);
            
            if sum(reasonable_idx) > 3
                integral_value = mean(pv_estimates(reasonable_idx));
                std_val = std(pv_estimates(reasonable_idx));
                fprintf('  奇点线去除法: %.6f ± %.6f (使用%d个估计值)\n', ...
                    integral_value, std_val, sum(reasonable_idx));
            else
                % 如果没有足够的合理估计，使用所有有效估计
                integral_value = mean(pv_estimates(valid_idx));
                std_val = std(pv_estimates(valid_idx));
                fprintf('  使用所有有效估计: %.6f ± %.6f\n', integral_value, std_val);
            end
        else
            fprintf('  警告: 无法获得收敛的主值估计\n');
            integral_value = 0;
        end
        
        % 方法2：如果方法1不稳定，尝试变换变量法
        % 利用 sin²(k₁/2) - sin²(k₂/2) = -sin((k₁+k₂)/2)·sin((k₁-k₂)/2)
        % 但这里不展开，使用方法1的结果
        
    % ========== 情况2：ω²LC ≠ 1（正常处理）==========
    else
        method_used{idx} = 'ω²LC≠1处理';
        
        % 分母: ω²LC·sin²(k₁/2) - sin²(k₂/2)
        denom = omega_sq_LC * sin(K1/2).^2 - sin(K2/2).^2;
        
        % 检查是否存在奇点（分母过零点）
        denom_min = min(denom(:));
        denom_max = max(denom(:));
        
        has_singularity = (denom_min * denom_max < 0);
        
        if has_singularity
            % 存在奇点，使用对称去除法
            
            % 奇点条件: |denom| < epsilon
            % 对于 ω²LC ≠ 1，奇点在
            % k₂ = ±2·arcsin(√(ω²LC)·sin(k₁/2)) 处
            
            eps_vals = logspace(-10, -3, 80);
            pv_estimates = zeros(size(eps_vals));
            
            for ei = 1:length(eps_vals)
                epsilon = eps_vals(ei);
                exclude_mask = abs(denom) < epsilon;
                
                % 控制排除区域大小
                if sum(exclude_mask(:)) > 0.3 * N_quad^2
                    pv_estimates(ei) = NaN;
                    continue;
                end
                
                % 去除小区域的积分
                valid_mask = ~exclude_mask;
                f_valid = num(valid_mask) ./ denom(valid_mask);
                pv_estimates(ei) = sum(f_valid(:)) * dA;
            end
            
            % 取收敛估计值
            valid_est = ~isnan(pv_estimates);
            if sum(valid_est) > 3
                % 选择后几个估计值的平均
                n_avg = min(10, sum(valid_est));
                last_valid = find(valid_est, n_avg, 'last');
                integral_value = mean(pv_estimates(last_valid));
            else
                % 如果所有epsilon都不合适，直接积分
                integral_value = sum(num(:) ./ denom(:)) * dA;
            end
        else
            % 无奇点，直接积分
            integral_value = sum(num(:) ./ denom(:)) * dA;
        end
    end
    
    % 存储结果
    integral_values(idx) = integral_value;
    ImZ_values(idx) = -omega_sq_LC / (4 * pi^2) * integral_value;
    
    % 显示进度
    if mod(idx, 200) == 0 || idx == 1 || idx == length(omega_sq_LC_vals) || abs(omega_sq_LC - 1.0) < 1e-15
        fprintf('  idx=%d, ω²LC=%.6f, Im[Z]=%.6f [%s]\n', ...
            idx, omega_sq_LC, ImZ_values(idx), method_used{idx});
    end
end

elapsed_time = toc;
fprintf('\n计算完成！总耗时: %.2f 秒\n', elapsed_time);

%% 输出 ω²LC = 1 的结果
idx_1 = find(abs(omega_sq_LC_vals - 1.0) < 1e-15, 1);
if ~isempty(idx_1)
    fprintf('\n========== ω²LC = 1 计算结果 ==========\n');
    fprintf('积分值 = %.8f\n', integral_values(idx_1));
    fprintf('Im[Z~(ω)] = %.8f\n', ImZ_values(idx_1));
end

%% ==================== 绘图 ====================

% 图1：完整范围
figure('Position', [100, 100, 1200, 500]);

% 子图1：积分值
subplot(1, 2, 1);
plot(omega_sq_LC_vals, integral_values, 'b-', 'LineWidth', 1.5);
xlabel('ω²LC', 'FontSize', 12);
ylabel('积分值（柯西主值）', 'FontSize', 12);
title(sprintf('柯西主值积分 (m=%d, n=%d)', m_fixed, n_fixed));
grid on;
hold on;
xline(1, 'r--', 'LineWidth', 1.5);
plot(omega_sq_LC_vals(idx_1), integral_values(idx_1), 'ro', ...
    'MarkerSize', 10, 'LineWidth', 2);
hold off;

% 子图2：Im[Z~(ω)]
subplot(1, 2, 2);
plot(omega_sq_LC_vals, ImZ_values, 'r-', 'LineWidth', 1.5);
xlabel('ω²LC', 'FontSize', 12);
ylabel('Im{Z~(ω)}', 'FontSize', 12);
title(sprintf('Im{{Z~(ω)}} (m=%d, n=%d)', m_fixed, n_fixed));
grid on;
hold on;
xline(1, 'r--', 'LineWidth', 1.5);
plot(omega_sq_LC_vals(idx_1), ImZ_values(idx_1), 'bo', ...
    'MarkerSize', 10, 'LineWidth', 2);
hold off;

% 图2：ω²LC=1附近的放大图
figure('Position', [200, 200, 800, 400]);
near_1 = abs(omega_sq_LC_vals - 1.0) < 0.0005;
plot(omega_sq_LC_vals(near_1), ImZ_values(near_1), 'b.-', ...
    'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('ω²LC', 'FontSize', 12);
ylabel('Im{Z~(ω)}', 'FontSize', 12);
title(sprintf('Im{{Z~(ω)}} 在ω²LC=1附近 (m=%d, n=%d)', m_fixed, n_fixed));
grid on;
hold on;
xline(1, 'r--', 'LineWidth', 1.5);
plot(1, ImZ_values(idx_1), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
text(1.00005, ImZ_values(idx_1), sprintf('Im[Z] = %.6f', ImZ_values(idx_1)), ...
    'FontSize', 10, 'BackgroundColor', 'w');
hold off;

%% 图3：不同m值的热力图
figure('Position', [300, 300, 800, 400]);

% 为多个m值计算（避免重复计算ω²LC=1的情况）
m_range = -5:1:5;
ImZ_matrix = zeros(length(m_range), length(omega_sq_LC_vals));

fprintf('\n计算多m值热力图...\n');
tic;

for mi = 1:length(m_range)
    m = m_range(mi);
    num_temp = 1 - cos(m*K1) .* cos(n_fixed*K2);
    
    for idx = 1:length(omega_sq_LC_vals)
        omega_sq_LC = omega_sq_LC_vals(idx);
        
        if abs(omega_sq_LC - 1.0) < 1e-15
            % ω²LC = 1的情况
            denom = sin(K1/2).^2 - sin(K2/2).^2;
            
            % 使用奇点线去除法
            delta = 0.05;  % 固定一个合适的delta
            near_main = abs(K1 - K2) < delta | abs(K1 - K2 - 2*pi) < delta | abs(K1 - K2 + 2*pi) < delta;
            near_anti = abs(K1 + K2) < delta | abs(K1 + K2 - 2*pi) < delta | abs(K1 + K2 + 2*pi) < delta;
            exclude_mask = near_main | near_anti;
            
            if sum(exclude_mask(:)) < 0.4 * N_quad^2
                valid_mask = ~exclude_mask;
                f_valid = num_temp(valid_mask) ./ denom(valid_mask);
                integral_value = sum(f_valid(:)) * dA;
            else
                integral_value = 0;
            end
        else
            % ω²LC ≠ 1的情况
            denom = omega_sq_LC * sin(K1/2).^2 - sin(K2/2).^2;
            
            if min(denom(:)) * max(denom(:)) < 0
                % 有奇点
                epsilon = 1e-4;  % 固定epsilon
                exclude_mask = abs(denom) < epsilon;
                if sum(exclude_mask(:)) < 0.3 * N_quad^2
                    valid_mask = ~exclude_mask;
                    f_valid = num_temp(valid_mask) ./ denom(valid_mask);
                    integral_value = sum(f_valid(:)) * dA;
                else
                    integral_value = sum(num_temp(:) ./ denom(:)) * dA;
                end
            else
                integral_value = sum(num_temp(:) ./ denom(:)) * dA;
            end
        end
        
        ImZ_matrix(mi, idx) = -omega_sq_LC / (4 * pi^2) * integral_value;
    end
    
    if mod(mi, 2) == 0
        fprintf('  m=%d 完成\n', m);
    end
end

fprintf('热力图计算完成，耗时: %.2f 秒\n', toc);

% 绘制
imagesc(omega_sq_LC_vals, m_range, ImZ_matrix);
colormap(jet(256));
clim_max = max(abs(ImZ_matrix(:)));
caxis([-clim_max, clim_max]);
colorbar;
xlabel('ω²LC', 'FontSize', 12);
ylabel('m', 'FontSize', 12);
title(sprintf('Im{{Z~(ω)}} 热力图 (n=%d)', n_fixed));
set(gca, 'YDir', 'normal');
hold on;
xline(1, 'k--', 'LineWidth', 2);
hold off;

%% 输出统计信息
fprintf('\n========== 统计信息 ==========\n');
fprintf('总点数: %d\n', length(omega_sq_LC_vals));
fprintf('ω²LC=1处的Im[Z]: %.6f\n', ImZ_values(idx_1));
fprintf('Im[Z]范围: [%.6f, %.6f]\n', min(ImZ_values), max(ImZ_values));
fprintf('Im[Z]均值: %.6f\n', mean(ImZ_values(~isnan(ImZ_values))));
fprintf('Im[Z]标准差: %.6f\n', std(ImZ_values(~isnan(ImZ_values))));

fprintf('\n========== 完成 ==========\n');