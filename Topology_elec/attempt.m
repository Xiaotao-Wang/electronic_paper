% 计算固定 n 时，Im[Z~(ω)] 随 m 和 ω²LC 的变化
% 使用高精度柯西主值积分 - 生成独立图形版本

clear; clc; close all;

%% 1. 参数设置
L = 1; C = 1;
n_fixed = 0;  % 固定 n (可以修改)

% 变量范围
omega_sq_LC_vals = linspace(0.999, 1.001, 1000);  % ω²LC 变化范围
m_vals = -5:0.1:5;  % m 变化范围

% 高精度积分设置
N_quad = 1000;  % 积分网格点数（平衡精度和速度）

fprintf('========== 高精度柯西主值积分计算 ==========\n');
fprintf('固定 n=%d，计算 Im[Z~(ω)] 随 (m, ω²LC) 的变化\n', n_fixed);
fprintf('积分网格: %d × %d (高斯-勒让德求积)\n', N_quad, N_quad);
fprintf('总计算点数: %d × %d = %d\n', length(m_vals), length(omega_sq_LC_vals), ...
    length(m_vals)*length(omega_sq_LC_vals));

% 生成高斯-勒让德节点和权重
[x, w] = lgwt(N_quad, -pi, pi);
[K1, K2] = meshgrid(x, x);
[W1, W2] = meshgrid(w, w);
weight_prod = W1 .* W2;

% 存储结果矩阵
ImZ_matrix = zeros(length(m_vals), length(omega_sq_LC_vals));

tic;
fprintf('\n开始计算...\n');

%% 对每个 (m, ω²LC) 组合计算柯西主值积分
for mi = 1:length(m_vals)
    m = m_vals(mi);
    
    % 计算分子（与 ω²LC 无关，但与 m 有关）
    num = 1 - cos(m*K1) .* cos(n_fixed*K2);
    
    for wi = 1:length(omega_sq_LC_vals)
        omega_sq_LC = omega_sq_LC_vals(wi);
        
        % 计算分母
        denom = omega_sq_LC * sin(K1/2).^2 - sin(K2/2).^2;
        
        % ========== 严格等于1时使用特殊处理 ==========
        if abs(omega_sq_LC - 1.0) < 1e-12
            % 使用小的虚部微扰
            eps_dist = 1e-9;
            f_matrix = num ./ (denom + 1i * eps_dist);
            f_real = real(f_matrix);
            
            % 处理奇异点区域
            denom_abs = abs(denom);
            singular_mask = denom_abs < 1e-4;
            
            if any(singular_mask(:))
                [sing_i, sing_j] = find(singular_mask);
                
                % 对每个奇异点进行局部平均处理
                for s = 1:min(length(sing_i), 500)
                    i0 = sing_i(s);
                    j0 = sing_j(s);
                    
                    radius = 10;
                    i_range = max(1, i0-radius):min(N_quad, i0+radius);
                    j_range = max(1, j0-radius):min(N_quad, j0+radius);
                    
                    if length(i_range) > 2 && length(j_range) > 2
                        denom_local = sin(K1(i_range, j_range)/2).^2 - sin(K2(i_range, j_range)/2).^2;
                        num_local = 1 - cos(m*K1(i_range, j_range)) .* cos(n_fixed*K2(i_range, j_range));
                        f_local = real(num_local ./ (denom_local + 1i * 1e-10));
                        
                        valid_mask = abs(denom_local) > 1e-5;
                        if sum(valid_mask(:)) > 10
                            f_local(~valid_mask) = mean(f_local(valid_mask));
                            f_real(i_range, j_range) = f_local;
                        end
                    end
                end
            end
            
            integral_value = sum(sum(f_real .* weight_prod));
            ImZ_matrix(mi, wi) = -omega_sq_LC / (4 * pi^2) * integral_value;
            continue;
        end
        
        % ========== 正常积分计算 ==========
        eps_dist = 1e-6;
        f_matrix = real(num ./ (denom + 1i * eps_dist));
        
        % 检测并处理奇异点
        denom_abs = abs(denom);
        singular_mask = denom_abs < 1e-3;
        
        if any(singular_mask(:))
            [sing_i, sing_j] = find(singular_mask);
            
            for s = 1:min(length(sing_i), 30)
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
        
        % 二维网格积分求和
        integral_value = sum(sum(f_matrix .* weight_prod));
        
        % 存储结果
        ImZ_matrix(mi, wi) = -omega_sq_LC / (4 * pi^2) * integral_value;
    end
    
    % 显示进度
    fprintf('进度: m=%d (%d/%d) 完成, 耗时: %.1f秒\n', m, mi, length(m_vals), toc);
end

elapsed_time = toc;
fprintf('\n计算完成！总耗时: %.2f 秒\n', elapsed_time);

%% 生成网格数据
[Omega_grid, M_grid] = meshgrid(omega_sq_LC_vals, m_vals);

%% ==================== 图1：三维Mesh图 ====================
figure('Name', '三维Mesh图', 'Position', [100, 100, 900, 700], 'Color', 'w');

% 创建三维曲面图
mesh(Omega_grid, M_grid, ImZ_matrix, 'FaceAlpha', 0.6);

% 美化图形
xlabel('$\omega^2 L C$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$m$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('$\mathrm{Im}[\tilde{Z}(\omega)]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('$\\mathrm{Im}[\\tilde{Z}(\\omega)]$ 三维分布 (固定 $n=%d$)', n_fixed), ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

% 颜色映射
colormap(jet(256));
colorbar;
caxis([-max(abs(ImZ_matrix(:))), max(abs(ImZ_matrix(:)))]);

% 视角设置
view(45, 30);
grid on;
box on;

% 添加光源效果
light('Position', [1, 1, 1], 'Style', 'infinite');
lighting gouraud;
material shiny;

% 保存图形
saveas(gcf, sprintf('3D_Mesh_n%d.png', n_fixed));
fprintf('图1已保存: 3D_Mesh_n%d.png\n', n_fixed);

%% ==================== 图2：m方向投影（固定omega） ====================
figure('Name', 'm方向投影', 'Position', [100, 100, 800, 500], 'Color', 'w');

% 选择几个代表性的omega值
omega_indices = [1, round(length(omega_sq_LC_vals)/3), round(length(omega_sq_LC_vals)/2), ...
                 round(2*length(omega_sq_LC_vals)/3), round(length(omega_sq_LC_vals)*0.9)];
omega_selected = omega_sq_LC_vals(omega_indices);

% 颜色映射
colors = lines(length(omega_indices));

hold on;
for i = 1:length(omega_indices)
    plot(m_vals, ImZ_matrix(:, omega_indices(i)), 'o-', ...
        'Color', colors(i,:), 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', sprintf('\\omega^2LC = %.3f', omega_selected(i)));
end
hold off;

xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$\mathrm{Im}[\tilde{Z}(\omega)]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('$\\mathrm{Im}[\\tilde{Z}(\\omega)]$ 随 $m$ 变化 (固定 $n=%d$)', n_fixed), ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

legend('Location', 'best', 'FontSize', 12);
grid on;
box on;

% 添加零线
hold on;
plot([min(m_vals), max(m_vals)], [0, 0], 'k--', 'LineWidth', 1.5);
hold off;

saveas(gcf, sprintf('Projection_m_n%d.png', n_fixed));
fprintf('图2已保存: Projection_m_n%d.png\n', n_fixed);

%% ==================== 图3：omega方向投影（固定m） ====================
figure('Name', 'omega方向投影', 'Position', [100, 100, 800, 500], 'Color', 'w');

% 选择几个代表性的m值
m_selected_indices = [1, round(length(m_vals)/4), round(length(m_vals)/2), ...
                      round(3*length(m_vals)/4), length(m_vals)];
m_selected = m_vals(m_selected_indices);

% 颜色映射
colors = lines(length(m_selected_indices));

hold on;
for i = 1:length(m_selected_indices)
    plot(omega_sq_LC_vals, ImZ_matrix(m_selected_indices(i), :), 'o-', ...
        'Color', colors(i,:), 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', sprintf('m = %.1f', m_selected(i)));
end
hold off;

xlabel('$\omega^2 L C$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$\mathrm{Im}[\tilde{Z}(\omega)]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('$\\mathrm{Im}[\\tilde{Z}(\\omega)]$ 随 $\\omega^2LC$ 变化 (固定 $n=%d$)', n_fixed), ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

legend('Location', 'best', 'FontSize', 12);
grid on;
box on;

% 标记 omega^2LC = 1 的位置
hold on;
yl = ylim;
plot([1, 1], yl, 'r--', 'LineWidth', 2, 'DisplayName', '\omega^2LC = 1');
hold off;

saveas(gcf, sprintf('Projection_omega_n%d.png', n_fixed));
fprintf('图3已保存: Projection_omega_n%d.png\n', n_fixed);

%% ==================== 图4：等高线图 ====================
figure('Name', '等高线图', 'Position', [100, 100, 800, 600], 'Color', 'w');

% 创建等高线图
contourf(Omega_grid, M_grid, ImZ_matrix, 50, 'LineStyle', 'none');
colormap(jet(256));
colorbar;

xlabel('$\omega^2 L C$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$m$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('$\\mathrm{Im}[\\tilde{Z}(\\omega)]$ 等高线图 (固定 $n=%d$)', n_fixed), ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

% 添加等值线
hold on;
contour(Omega_grid, M_grid, ImZ_matrix, 10, 'k-', 'LineWidth', 1);
% 标记零线
contour(Omega_grid, M_grid, ImZ_matrix, [0, 0], 'k-', 'LineWidth', 2);
% 标记 omega^2LC = 1
plot([1, 1], [min(m_vals), max(m_vals)], 'r--', 'LineWidth', 2);
hold off;

grid on;
box on;

saveas(gcf, sprintf('Contour_n%d.png', n_fixed));
fprintf('图4已保存: Contour_n%d.png\n', n_fixed);

%% ==================== 图5：瀑布图（所有m值的omega方向投影） ====================
figure('Name', '瀑布图', 'Position', [100, 100, 1000, 700], 'Color', 'w');

% 创建瀑布图
waterfall(Omega_grid, M_grid, ImZ_matrix);

xlabel('$\omega^2 L C$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$m$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('$\mathrm{Im}[\tilde{Z}(\omega)]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('$\\mathrm{Im}[\\tilde{Z}(\\omega)]$ 瀑布图 (固定 $n=%d$)', n_fixed), ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

colormap(jet(256));
colorbar;
view(-30, 45);
grid on;
box on;

saveas(gcf, sprintf('Waterfall_n%d.png', n_fixed));
fprintf('图5已保存: Waterfall_n%d.png\n', n_fixed);

%% ==================== 图6：特定omega^2LC=1处，所有m的ImZ ====================
figure('Name', 'omega²LC=1 处截面', 'Position', [100, 100, 800, 500], 'Color', 'w');

% 找到最接近 omega^2LC=1 的索引
[~, idx_1] = min(abs(omega_sq_LC_vals - 1.0));
ImZ_at_one = ImZ_matrix(:, idx_1);

% 绘制
plot(m_vals, ImZ_at_one, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$\mathrm{Im}[\tilde{Z}(\omega)]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('$\\mathrm{Im}[\\tilde{Z}(\\omega)]$ 在 $\\omega^2LC=1$ 处随 $m$ 变化 (固定 $n=%d$)', n_fixed), ...
    'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');

grid on;
box on;

% 标记零点
hold on;
plot([min(m_vals), max(m_vals)], [0, 0], 'k--', 'LineWidth', 1.5);
hold off;

% 标记特殊点
[max_val, max_idx] = max(ImZ_at_one);
[min_val, min_idx] = min(ImZ_at_one);
hold on;
plot(m_vals(max_idx), max_val, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
plot(m_vals(min_idx), min_val, 'g*', 'MarkerSize', 15, 'LineWidth', 2);
text(m_vals(max_idx)+0.2, max_val, sprintf('Max: %.3f', max_val), 'FontSize', 10, 'Color', 'r');
text(m_vals(min_idx)+0.2, min_val, sprintf('Min: %.3f', min_val), 'FontSize', 10, 'Color', 'g');
hold off;

saveas(gcf, sprintf('At_omega2LC1_n%d.png', n_fixed));
fprintf('图6已保存: At_omega2LC1_n%d.png\n', n_fixed);

%% 输出统计信息
fprintf('\n========== 统计信息 ==========\n');
fprintf('固定 n = %d\n', n_fixed);
fprintf('m 范围: %.1f 到 %.1f (共 %d 个点)\n', min(m_vals), max(m_vals), length(m_vals));
fprintf('ω²LC 范围: %.3f 到 %.3f (共 %d 个点)\n', min(omega_sq_LC_vals), max(omega_sq_LC_vals), length(omega_sq_LC_vals));
fprintf('积分网格: %d × %d\n', N_quad, N_quad);
fprintf('Im[Z] 全局范围: %.4f 到 %.4f\n', min(ImZ_matrix(:)), max(ImZ_matrix(:)));
fprintf('Im[Z] 均值: %.4f\n', mean(ImZ_matrix(:)));
fprintf('Im[Z] 标准差: %.4f\n', std(ImZ_matrix(:)));
fprintf('在 ω²LC=1 处 Im[Z] 范围: %.4f 到 %.4f\n', min(ImZ_at_one), max(ImZ_at_one));

%% 自定义高斯-勒让德求积节点和权重函数
function [x, w] = lgwt(N, a, b)
    % 初始猜测（使用切比雪夫节点）
    x = cos(pi * (4*(1:N)' - 1) / (4*N + 2));
    epsilon = 1e-15;
    
    % Newton-Raphson迭代
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
    
    % 计算权重
    w = 2 ./ ((1 - x.^2) .* P_der_curr.^2);
    
    % 映射到区间 [a, b]
    x = 0.5 * (a + b) + 0.5 * (b - a) * x;
    w = 0.5 * (b - a) * w;
end