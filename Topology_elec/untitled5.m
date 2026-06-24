% 使用变量替换计算广义积分
% I(m,n) = ∫∫ [1-cos(mk1+nk2)] / [alpha*sin²(k1/2)-sin²(k2/2)] dk1 dk2
% 令 u = k1/2, v = k2/2

clear; clc; close all;

%% 1. 参数设置
alpha = 2;  % 可以修改

% 范围
m_vals = -10:1:10;
n_vals = -10:1:10;

fprintf('========== 使用变量替换计算积分 ==========\n');
fprintf('alpha = %.2f\n', alpha);
fprintf('令 u = k1/2, v = k2/2\n');
fprintf('I(m,n) = 4 * ∫∫ [1-cos(2mu+2nv)] / [alpha*sin²(u)-sin²(v)] du dv\n\n');

% 积分网格
N_quad = 400;
u = linspace(-pi/2+1e-10, pi/2-1e-10, N_quad);
v = linspace(-pi/2+1e-10, pi/2-1e-10, N_quad);
[du, dv] = meshgrid(u, v);

% 分母: alpha*sin²(u) - sin²(v)
denom = alpha * sin(du).^2 - sin(dv).^2;

% 奇点检测
threshold = 1e-12;
singular_idx = abs(denom) < threshold;

% 存储结果
I_matrix = zeros(length(m_vals), length(n_vals));

fprintf('网格点数: %d × %d\n', N_quad, N_quad);
fprintf('总计算量: %d 个积分点\n\n', length(m_vals)*length(n_vals));

total = length(m_vals) * length(n_vals);
progress = 0;
tic;

for mi = 1:length(m_vals)
    m = m_vals(mi);
    for ni = 1:length(n_vals)
        n = n_vals(ni);
        
        progress = progress + 1;
        if mod(progress, 50) == 0
            fprintf('进度: %d/%d (%.1f%%)\n', progress, total, progress/total*100);
        end
        
        if m == 0 && n == 0
            I_matrix(mi, ni) = 0;
            continue;
        end
        
        % 计算被积函数 (变量替换后多一个因子4)
        theta = 2 * (m * du + n * dv);
        numerator = 1 - cos(theta);
        integrand = numerator ./ denom;
        
        % 处理奇点
        if any(singular_idx(:))
            se = strel('disk', 3);
            singular_idx_expanded = imdilate(singular_idx, se);
            non_singular = ~singular_idx_expanded;
            
            if sum(non_singular(:)) > 0
                avg_val = mean(integrand(non_singular));
                integrand(singular_idx_expanded) = avg_val;
            else
                integrand(singular_idx_expanded) = 0;
            end
        end
        
        % 二重积分 (乘以雅可比行列式4)
        I_matrix(mi, ni) = 4 * trapz(v, trapz(u, integrand, 1));
    end
end

elapsed_time = toc;
fprintf('\n计算完成！总耗时: %.2f 秒\n', elapsed_time);

%% 2. 输出结果
fprintf('\n========== 积分值 (alpha = %.2f) ==========\n', alpha);

% 找出非零值
[m_idx, n_idx] = find(abs(I_matrix) > 1e-4);
non_zero_count = length(m_idx);

fprintf('非零值总数: %d / %d\n', non_zero_count, length(m_vals)*length(n_vals));
fprintf('稀疏度: %.2f%%\n', non_zero_count/(length(m_vals)*length(n_vals))*100);

% 显示小范围矩阵 (m,n = -5 到 5)
fprintf('\n小范围矩阵 (m,n = -5 到 5):\n');
fprintf('m\\n\t');
for n = -5:5
    fprintf('%6d\t', n);
end
fprintf('\n');

for m = -5:5
    fprintf('%d\t', m);
    for n = -5:5
        idx_m = find(m_vals == m);
        idx_n = find(n_vals == n);
        if ~isempty(idx_m) && ~isempty(idx_n)
            fprintf('%6.2f\t', I_matrix(idx_m, idx_n));
        end
    end
    fprintf('\n');
end

%% 3. 绘制结果
figure('Position', [50, 50, 1200, 800]);

% 热力图
subplot(2, 3, 1);
display_mat = I_matrix;
display_mat(abs(display_mat) < 1e-4) = 0;
imagesc(n_vals, m_vals, display_mat);
colorbar;
xlabel('n', 'FontSize', 12);
ylabel('m', 'FontSize', 12);
title(sprintf('I(m,n) 热力图 (alpha = %.2f)', alpha), 'FontSize', 13);
set(gca, 'YDir', 'normal');
axis square;

% 3D图
subplot(2, 3, 2);
[M, N] = meshgrid(m_vals, n_vals);
display_3d = I_matrix;
display_3d(abs(display_3d) < 1e-4) = NaN;
surf(M, N, display_3d', 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xlabel('m', 'FontSize', 11);
ylabel('n', 'FontSize', 11);
zlabel('I(m,n)', 'FontSize', 11);
title('3D 结构', 'FontSize', 13);
colormap('jet');
colorbar;
grid on;
view(45, 30);

% 沿m方向 (n=0)
subplot(2, 3, 3);
idx_n0 = find(n_vals == 0);
plot(m_vals, I_matrix(idx_n0, :), 'bo-', 'LineWidth', 2);
xlabel('m', 'FontSize', 12);
ylabel('I(m,0)', 'FontSize', 12);
title('沿 m 方向 (n=0)', 'FontSize', 13);
grid on;
xlim([-10, 10]);

% 沿n方向 (m=0)
subplot(2, 3, 4);
idx_m0 = find(m_vals == 0);
plot(n_vals, I_matrix(idx_m0, :), 'ro-', 'LineWidth', 2);
xlabel('n', 'FontSize', 12);
ylabel('I(0,n)', 'FontSize', 12);
title('沿 n 方向 (m=0)', 'FontSize', 13);
grid on;
xlim([-10, 10]);

% 统计
subplot(2, 3, 5);
text(0.1, 0.9, sprintf('alpha = %.2f', alpha), 'FontSize', 13, 'FontWeight', 'bold');
text(0.1, 0.7, sprintf('非零值: %d', non_zero_count), 'FontSize', 11);
text(0.1, 0.5, sprintf('稀疏度: %.2f%%', non_zero_count/(length(m_vals)*length(n_vals))*100), 'FontSize', 11);
if non_zero_count > 0
    text(0.1, 0.3, sprintf('最大值: %.4f', max(I_matrix(:))), 'FontSize', 11);
    text(0.1, 0.1, sprintf('最小值: %.4f', min(I_matrix(:))), 'FontSize', 11);
end
axis off;
title('统计信息', 'FontSize', 13);

% 非零值分布
subplot(2, 3, 6);
if non_zero_count > 0
    M_nonzero = m_vals(m_idx);
    N_nonzero = n_vals(n_idx);
    Z_nonzero = I_matrix(sub2ind(size(I_matrix), m_idx, n_idx));
    
    scatter3(M_nonzero, N_nonzero, Z_nonzero, 40, Z_nonzero, 'filled');
    xlabel('m', 'FontSize', 10);
    ylabel('n', 'FontSize', 10);
    zlabel('I(m,n)', 'FontSize', 10);
    title('非零值分布', 'FontSize', 12);
    colormap('jet');
    colorbar;
    grid on;
    view(45, 30);
end

%% 4. 保存
save(sprintf('Integral_alpha_%.2f_new.mat', alpha), 'm_vals', 'n_vals', 'I_matrix', 'alpha');
fprintf('\n数据已保存\n');