% 计算阻抗函数 Z~(ω) 的虚部（使用柯西主值积分）
% 被积函数分母: omega^2*L*C*sin(k1/2)^2 - sin(k2/2)^2

clear; clc; close all;

%% 1. 参数设置
L = 1; C = 1; omega = 0.1;
omega_sq_LC = omega^2 * L * C;  % omega²LC

% 扩大范围
m_vals = -8:1:8;
n_vals = -8:1:8;

% 高精度积分设置
N_quad = 1000;  % 积分网格点数
tol = 1e-10;    % 奇点检测阈值

fprintf('开始计算 Im[Z~(ω)] 数值积分（柯西主值）...\n');
fprintf('分母 = ω²LC·sin²(k₁/2) - sin²(k₂/2)\n');
fprintf('ω²LC = %.4f\n', omega_sq_LC);
fprintf('网格点数: %d × %d\n', N_quad, N_quad);
fprintf('m 范围: %d 到 %d, n 范围: %d 到 %d\n', ...
    min(m_vals), max(m_vals), min(n_vals), max(n_vals));
fprintf('总计算量: %d 个积分点\n\n', length(m_vals)*length(n_vals));

% 存储结果
ImZ_matrix = zeros(length(m_vals), length(n_vals));
total_points = length(m_vals) * length(n_vals);
progress = 0;

tic;

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
        
        % 使用柯西主值积分
        integral_value = cauchy_principal_value_2d(m, n, N_quad, omega_sq_LC, tol);
        
        % Z~(ω) 的虚部
        ImZ_matrix(mi, ni) = -omega * L / (4 * pi^2) * integral_value;
    end
end

elapsed_time = toc;
fprintf('\n计算完成！耗时: %.2f 秒\n', elapsed_time);

%% 被积函数
function f = integrand_func(k1, k2, m, n, omega_sq_LC)
    % 被积函数: [1 - cos(m*k1 + n*k2)] / [omega²LC*sin²(k1/2) - sin²(k2/2)]
    theta = m * k1 + n * k2;
    numerator = 1 - cos(theta);
    denom = omega_sq_LC * sin(k1/2).^2 - sin(k2/2).^2;
    
    % 防止除零
    if abs(denom) < 1e-15
        f = 0;
    else
        f = numerator ./ denom;
    end
end

%% 改进的柯西主值积分函数（二维）
function pv = cauchy_principal_value_2d(m, n, N, omega_sq_LC, tol)
    % 计算二维柯西主值积分
    % PV ∫∫ [1 - cos(m*k1 + n*k2)] / [omega²LC*sin²(k1/2) - sin²(k2/2)] dk1 dk2
    
    % 积分区间 [-π, π]
    a = -pi; b = pi;
    
    % 使用高斯-勒让德积分提高精度
    [x, w] = lgwt(N, a, b);
    
    % 预计算 sin² 值以提高效率
    sin2_k1 = sin(x/2).^2;
    sin2_k2 = sin(x/2).^2;
    
    % 初始化积分值
    pv = 0;
    
    % 使用二维积分，处理奇点
    for i = 1:N
        k1 = x(i);
        s1 = sin2_k1(i);
        
        for j = 1:N
            k2 = x(j);
            
            % 计算分母
            denom = omega_sq_LC * s1 - sin2_k2(j);
            
            % 检查是否为奇点
            if abs(denom) < tol
                % 对于奇点，使用对称邻域平均
                % 这里采用更精细的处理：在奇点附近进行自适应积分
                pv = pv + handle_singularity(k1, k2, m, n, omega_sq_LC, tol, w(i), w(j));
            else
                % 正常点
                theta = m * k1 + n * k2;
                numerator = 1 - cos(theta);
                f = numerator / denom;
                pv = pv + w(i) * w(j) * f;
            end
        end
    end
end

%% 处理奇点的子函数
function val = handle_singularity(k1, k2, m, n, omega_sq_LC, tol, wi, wj)
    % 在奇点附近使用小邻域积分
    % 奇点满足: omega²LC*sin²(k1/2) - sin²(k2/2) = 0
    
    % 使用小邻域半径
    delta = 1e-4;
    n_sub = 20;  % 子网格点数
    
    % 在奇点周围的小区域内进行数值积分
    val = 0;
    
    % 使用极坐标在奇点周围积分
    for r = 1:n_sub
        for theta_idx = 1:n_sub
            % 极坐标参数
            rho = delta * r / n_sub;
            theta_angle = 2 * pi * (theta_idx - 0.5) / n_sub;
            
            % 计算偏移量
            dk1 = rho * cos(theta_angle);
            dk2 = rho * sin(theta_angle);
            
            % 采样点
            k1_s = k1 + dk1;
            k2_s = k2 + dk2;
            
            % 检查是否在积分范围内
            if abs(k1_s) <= pi && abs(k2_s) <= pi
                % 计算被积函数
                denom = omega_sq_LC * sin(k1_s/2)^2 - sin(k2_s/2)^2;
                theta = m * k1_s + n * k2_s;
                numerator = 1 - cos(theta);
                
                if abs(denom) > tol
                    f = numerator / denom;
                    % 极坐标积分权重
                    val = val + f * rho * delta/n_sub * 2*pi/n_sub;
                end
            end
        end
    end
    
    % 乘以高斯权重
    val = val * wi * wj;
end

%% 高斯-勒让德求积节点和权重
function [x, w] = lgwt(N, a, b)
    % 高斯-勒让德求积公式
    % N: 节点数, [a,b]: 积分区间
    % 返回: x - 节点, w - 权重
    
    % 初始猜测
    x = cos(pi * (4*(1:N)' - 1) / (4*N + 2));
    
    % 牛顿迭代
    epsilon = 1e-15;
    for iter = 1:100
        P = zeros(N, 1);
        P_der = zeros(N, 1);
        
        % 勒让德多项式递推
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
        
        P = P_curr;
        P_der = P_der_curr;
        
        % 更新
        dx = P ./ P_der;
        x = x - dx;
        
        if max(abs(dx)) < epsilon
            break;
        end
    end
    
    % 计算权重
    w = 2 ./ ((1 - x.^2) .* P_der.^2);
    
    % 映射到 [a, b]
    x = 0.5 * (a + b) + 0.5 * (b - a) * x;
    w = 0.5 * (b - a) * w;
end

%% 显示结果
fprintf('\n========== 非零值检测 ==========\n');
fprintf('阈值: 1e-6\n\n');

% 找出所有非零值
[m_idx, n_idx] = find(abs(ImZ_matrix) > 1e-6);
non_zero_count = length(m_idx);

if non_zero_count == 0
    fprintf('未检测到非零值\n');
else
    fprintf('检测到 %d 个非零值:\n', non_zero_count);
    fprintf('m\tn\tIm[Z~(ω)]\n');
    fprintf('------------------------\n');
    
    % 存储非零值用于分析
    non_zero_values = [];
    for i = 1:non_zero_count
        m_val = m_vals(m_idx(i));
        n_val = n_vals(n_idx(i));
        z_val = ImZ_matrix(m_idx(i), n_idx(i));
        fprintf('%3d\t%3d\t%.8f\n', m_val, n_val, z_val);
        non_zero_values = [non_zero_values; m_val, n_val, z_val];
    end
    
    % 分析模式
    fprintf('\n========== 模式分析 ==========\n');
    fprintf('按 |m|, |n| 分组:\n');
    fprintf('|m|\t|n|\t值\n');
    fprintf('------------------------\n');
    
    for i = 1:non_zero_count
        m_abs = abs(non_zero_values(i, 1));
        n_abs = abs(non_zero_values(i, 2));
        z_val = non_zero_values(i, 3);
        fprintf('%d\t%d\t%.8f\n', m_abs, n_abs, z_val);
    end
end

%% 绘制沿 m 方向的子图（固定 n）
figure('Name', '柯西主值积分结果', 'Position', [100, 100, 1200, 700]);

% 选择关键 n 值
n_fixed = [-3, -2, -1, 0, 1, 2, 3];
colors = lines(length(n_fixed));

hold on;
for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    idx_n = find(n_vals == n);
    if ~isempty(idx_n)
        ImZ_values = ImZ_matrix(:, idx_n)';
        
        plot(m_vals, ImZ_values, 'o-', 'LineWidth', 2, 'Color', colors(ni,:), ...
            'MarkerSize', 6, 'DisplayName', sprintf('n = %d', n));
    end
end

xlabel('m', 'FontSize', 13);
ylabel('Im{Z~(ω)}', 'FontSize', 13);
title(sprintf('柯西主值积分: 沿 m 方向变化（固定 n），ω²LC = %.2f', omega_sq_LC), 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
grid on;
xlim([min(m_vals), max(m_vals)]);
set(gca, 'FontSize', 11);
hold off;

%% 三维曲面图
figure('Name', '三维视图', 'Position', [150, 150, 1000, 700]);

[M, N] = meshgrid(m_vals, n_vals);

% 使用散点图显示所有非零值
nonzero_idx = abs(ImZ_matrix) > 1e-6;
M_nonzero = M(nonzero_idx);
N_nonzero = N(nonzero_idx);
Z_nonzero = ImZ_matrix(nonzero_idx);

if ~isempty(Z_nonzero)
    scatter3(M_nonzero, N_nonzero, Z_nonzero, 150, Z_nonzero, 'filled');
    xlabel('m', 'FontSize', 12);
    ylabel('n', 'FontSize', 12);
    zlabel('Im{Z~(ω)}', 'FontSize', 12);
    title(sprintf('柯西主值积分: Im{Z~(ω)} 非零值分布 (共 %d 个点)', length(Z_nonzero)), 'FontSize', 13);
    colormap(jet);
    colorbar;
    grid on;
    view(45, 30);
    
    % 添加标注
    for i = 1:length(M_nonzero)
        text(M_nonzero(i), N_nonzero(i), Z_nonzero(i), ...
            sprintf('(%d,%d)', M_nonzero(i), N_nonzero(i)), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
            'FontSize', 9);
    end
else
    text(0, 0, 0, '无非零值', 'FontSize', 14, 'HorizontalAlignment', 'center');
    xlabel('m', 'FontSize', 12);
    ylabel('n', 'FontSize', 12);
    zlabel('Im{Z~(ω)}', 'FontSize', 12);
    title('柯西主值积分: Im{Z~(ω)} 三维视图', 'FontSize', 13);
    grid on;
end

%% 热力图
figure('Name', '热力图', 'Position', [200, 200, 800, 700]);

imagesc(n_vals, m_vals, ImZ_matrix);
colormap(redbluecmap);
colorbar;
xlabel('n', 'FontSize', 12);
ylabel('m', 'FontSize', 12);
title(sprintf('柯西主值积分: Im{Z~(ω)} 热力图, ω²LC = %.2f', omega_sq_LC), 'FontSize', 13);
set(gca, 'YDir', 'normal');
axis square;

% 标记非零值
hold on;
for mi = 1:length(m_vals)
    for ni = 1:length(n_vals)
        if abs(ImZ_matrix(mi, ni)) > 1e-6
            val = ImZ_matrix(mi, ni);
            if abs(val) > 0.5
                text_color = 'w';
            else
                text_color = 'k';
            end
            text(n_vals(ni), m_vals(mi), sprintf('%.2f', val), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', text_color, 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
end
hold off;

%% 自定义红蓝配色函数
function cmap = redbluecmap
    n = 256;
    cmap = zeros(n, 3);
    half = floor(n/2);
    
    % 蓝色到白色（负值区域）
    for i = 1:half
        t = (i-1)/(half-1);
        cmap(i, :) = [t, t, 1];
    end
    
    % 白色到红色（正值区域）
    for i = half+1:n
        t = (i-half-1)/(n-half-1);
        cmap(i, :) = [1, 1-t, 1-t];
    end
end

%% 保存数据
save('Impedance_Cauchy_Principal_modified.mat', 'm_vals', 'n_vals', 'ImZ_matrix', 'omega', 'L', 'C', 'omega_sq_LC');
fprintf('\n数据已保存到 Impedance_Cauchy_Principal_modified.mat\n');

%% 输出统计信息
fprintf('\n========== 统计信息 ==========\n');
fprintf('总点数: %d\n', length(m_vals)*length(n_vals));
fprintf('非零值数量: %d\n', non_zero_count);
if non_zero_count > 0
    fprintf('稀疏度: %.2f%%\n', non_zero_count/(length(m_vals)*length(n_vals))*100);
    fprintf('最大值: %.8f\n', max(abs(ImZ_matrix(:))));
    fprintf('最小值（非零）: %.8f\n', min(abs(ImZ_matrix(abs(ImZ_matrix)>1e-6))));
end