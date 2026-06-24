% 计算阻抗函数 Z~(ω) 的虚部（扩大范围，使用柯西主值积分）
% 使用柯西主值积分处理奇点

clear; clc; close all;

%% 1. 参数设置
L = 1; C = 1; omega = 1;

% 扩大范围
m_vals = -8:1:8;          % m 从 -8 到 8
n_vals = -8:1:8;          % n 从 -8 到 8

% 高精度积分网格
N_quad = 800;             % 提高精度
fprintf('开始计算 Im[Z~(ω)] 数值积分（柯西主值）...\n');
fprintf('ω²LC = 1, 网格点数: %d × %d\n', N_quad, N_quad);
fprintf('m 范围: %d 到 %d, n 范围: %d 到 %d\n', ...
    min(m_vals), max(m_vals), min(n_vals), max(n_vals));
fprintf('总计算量: %d 个积分点\n\n', length(m_vals)*length(n_vals));

% 存储结果
ImZ_matrix = zeros(length(m_vals), length(n_vals));

% 计数器
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
        
        % m=n=0 时特殊处理
        if m == 0 && n == 0
            ImZ_matrix(mi, ni) = 0;
            continue;
        end
        
        % === 使用柯西主值积分 ===
        integral_value = cauchy_principal_value_2d(m, n, N_quad);
        
        % Z~(ω) 的虚部
        ImZ_matrix(mi, ni) = -omega * L / (4 * pi^2) * integral_value;
    end
end

elapsed_time = toc;
fprintf('\n计算完成！耗时: %.2f 秒\n', elapsed_time);

%% 被积函数
function f = integrand_func(k1, k2, m, n)
    % 被积函数: [1 - cos(m*k1 + n*k2)] / [sin(k1/2)^2 - sin(k2/2)^2]
    theta = m * k1 + n * k2;
    numerator = 1 - cos(theta);
    denom = sin(k1/2).^2 - sin(k2/2).^2;
    
    % 防止除零（在柯西主值中，奇点处用极限值）
    if abs(denom) < 1e-15
        % 在奇点附近，使用洛必达法则或对称极限
        % 对于奇点，分子也为0，极限值取决于方向
        f = 0;
    else
        f = numerator ./ denom;
    end
end

%% 柯西主值积分函数（二维）- 改进版
function pv = cauchy_principal_value_2d(m, n, N)
    % 计算二维柯西主值积分
    % PV ∫∫ [1 - cos(m*k1 + n*k2)] / [sin(k1/2)^2 - sin(k2/2)^2] dk1 dk2
    
    % 积分区间
    a = -pi; b = pi;
    
    % 使用奇数个点保证对称性
    N_odd = N;
    if mod(N_odd, 2) == 0
        N_odd = N_odd + 1;
    end
    
    % 生成网格
    k1 = linspace(a, b, N_odd);
    k2 = linspace(a, b, N_odd);
    [K1, K2] = meshgrid(k1, k2);
    
    % 计算被积函数
    theta = m * K1 + n * K2;
    numerator = 1 - cos(theta);
    denom = sin(K1/2).^2 - sin(K2/2).^2;
    
    % 检测奇点（分母为零的点）
    singular_threshold = 1e-10;
    singular_mask = abs(denom) < singular_threshold;
    
    % 初始化被积函数矩阵
    F = zeros(size(K1));
    
    % === 方法1: 使用对称挖去法 ===
    % 对所有非奇点正常计算
    F(~singular_mask) = numerator(~singular_mask) ./ denom(~singular_mask);
    
    % 对奇点使用对称邻域平均
    if any(singular_mask(:))
        % 对每个奇点，使用其对称方向的邻域值
        for i = 1:N_odd
            for j = 1:N_odd
                if singular_mask(i, j)
                    % 收集周围非奇点的值
                    neighbor_vals = [];
                    
                    % 检查8个方向
                    for di = -2:2
                        for dj = -2:2
                            if di == 0 && dj == 0
                                continue;
                            end
                            i2 = i + di;
                            j2 = j + dj;
                            if i2 >= 1 && i2 <= N_odd && j2 >= 1 && j2 <= N_odd
                                if ~singular_mask(i2, j2)
                                    % 计算该点的被积函数值
                                    theta_nb = m * K1(i2, j2) + n * K2(i2, j2);
                                    num_nb = 1 - cos(theta_nb);
                                    denom_nb = sin(K1(i2, j2)/2).^2 - sin(K2(i2, j2)/2).^2;
                                    if abs(denom_nb) > singular_threshold
                                        neighbor_vals = [neighbor_vals, num_nb / denom_nb];
                                    end
                                end
                            end
                        end
                    end
                    
                    if ~isempty(neighbor_vals)
                        % 使用邻域的平均值
                        F(i, j) = mean(neighbor_vals);
                    else
                        % 如果周围都是奇点，使用更远的点
                        for radius = 3:5
                            neighbor_vals = [];
                            for di = -radius:radius
                                for dj = -radius:radius
                                    if di == 0 && dj == 0
                                        continue;
                                    end
                                    i2 = i + di;
                                    j2 = j + dj;
                                    if i2 >= 1 && i2 <= N_odd && j2 >= 1 && j2 <= N_odd
                                        if ~singular_mask(i2, j2)
                                            theta_nb = m * K1(i2, j2) + n * K2(i2, j2);
                                            num_nb = 1 - cos(theta_nb);
                                            denom_nb = sin(K1(i2, j2)/2).^2 - sin(K2(i2, j2)/2).^2;
                                            if abs(denom_nb) > singular_threshold
                                                neighbor_vals = [neighbor_vals, num_nb / denom_nb];
                                            end
                                        end
                                    end
                                end
                            end
                            if ~isempty(neighbor_vals)
                                F(i, j) = mean(neighbor_vals);
                                break;
                            end
                        end
                    end
                    
                    % 如果仍然没有邻域值，设置为0
                    if F(i, j) == 0 && singular_mask(i, j)
                        F(i, j) = 0;
                    end
                end
            end
        end
    end
    
    % 使用梯形法则进行二重积分
    pv = trapz(k2, trapz(k1, F, 1));
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
title(sprintf('柯西主值积分: 沿 m 方向变化（固定 n），ω²LC = 1'), 'FontSize', 14);
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
    colormap(jet);  % 三维图使用jet
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

%% 热力图 - 修改为红色系
figure('Name', '热力图', 'Position', [200, 200, 800, 700]);

imagesc(n_vals, m_vals, ImZ_matrix);

% === 修改热力图的颜色映射：将黄色改为红色 ===
% 方法1：使用 red-blue colormap（推荐）
colormap(redbluecmap);  % 自定义红蓝配色

% 或者方法2：使用内置的 redbluecmap（如果可用）
% colormap(redblue);  % MATLAB R2019b+

% 或者方法3：自定义红色系（从白到红）
% n_colors = 256;
% red_cmap = [linspace(1, 1, n_colors)', linspace(1, 0, n_colors)', linspace(1, 0, n_colors)'];
% colormap(red_cmap);

% 或者方法4：使用 hot colormap（从黑到红到黄到白）
% colormap(hot);

% 或者方法5：使用 parula 但调整
% colormap(parula);

colorbar;
xlabel('n', 'FontSize', 12);
ylabel('m', 'FontSize', 12);
title('柯西主值积分: Im{Z~(ω)} 热力图 (红蓝配色)', 'FontSize', 13);
set(gca, 'YDir', 'normal');
axis square;

% 标记非零值
hold on;
for mi = 1:length(m_vals)
    for ni = 1:length(n_vals)
        if abs(ImZ_matrix(mi, ni)) > 1e-6
            % 根据值的大小决定文字颜色（白色或黑色）
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
    % 创建红蓝配色：红色表示正值，蓝色表示负值
    n = 256;
    % 从蓝色到白色到红色
    cmap = zeros(n, 3);
    half = floor(n/2);
    
    % 蓝色到白色（负值区域）
    for i = 1:half
        t = (i-1)/(half-1);
        cmap(i, :) = [t, t, 1];  % 从 [0,0,1] 到 [1,1,1]
    end
    
    % 白色到红色（正值区域）
    for i = half+1:n
        t = (i-half-1)/(n-half-1);
        cmap(i, :) = [1, 1-t, 1-t];  % 从 [1,1,1] 到 [1,0,0]
    end
end

%% 保存数据
save('Impedance_Cauchy_Principal.mat', 'm_vals', 'n_vals', 'ImZ_matrix', 'omega', 'L', 'C');
fprintf('\n数据已保存到 Impedance_Cauchy_Principal.mat\n');

%% 输出统计信息
fprintf('\n========== 统计信息 ==========\n');
fprintf('总点数: %d\n', length(m_vals)*length(n_vals));
fprintf('非零值数量: %d\n', non_zero_count);
if non_zero_count > 0
    fprintf('稀疏度: %.2f%%\n', non_zero_count/(length(m_vals)*length(n_vals))*100);
    fprintf('最大值: %.8f\n', max(abs(ImZ_matrix(:))));
    fprintf('最小值（非零）: %.8f\n', min(abs(ImZ_matrix(abs(ImZ_matrix)>1e-6))));
end