% 直接数值计算原始积分 (扩大范围)
% I(m,n) = ∫∫ [1-cos(mk1+nk2)] / [sin²(k1/2)-sin²(k2/2)] dk1 dk2

clear; clc; close all;

%% 1. 参数设置
% 积分范围
m_vals = -20:1:20;
n_vals = -20:1:20;

% 高精度积分网格
N_quad = 800;
fprintf('========== 直接数值计算原始积分 ==========\n');
fprintf('网格点数: %d × %d\n', N_quad, N_quad);
fprintf('m 范围: %d 到 %d (共 %d 个)\n', min(m_vals), max(m_vals), length(m_vals));
fprintf('n 范围: %d 到 %d (共 %d 个)\n', min(n_vals), max(n_vals), length(n_vals));
fprintf('总计算量: %d 个积分点\n\n', length(m_vals)*length(n_vals));

% 积分网格
k1 = linspace(-pi+1e-10, pi-1e-10, N_quad);
k2 = linspace(-pi+1e-10, pi-1e-10, N_quad);
[dk1, dk2] = meshgrid(k1, k2);

% 预计算分母
denom = sin(dk1/2).^2 - sin(dk2/2).^2;
threshold = 1e-12;
singular_idx = abs(denom) < threshold;

% 存储结果
I_matrix = zeros(length(m_vals), length(n_vals));

% 进度条
total = length(m_vals) * length(n_vals);
progress = 0;
tic;

for mi = 1:length(m_vals)
    m = m_vals(mi);
    for ni = 1:length(n_vals)
        n = n_vals(ni);
        
        progress = progress + 1;
        if mod(progress, 100) == 0
            fprintf('进度: %d/%d (%.1f%%)  耗时: %.1f秒\n', ...
                progress, total, progress/total*100, toc);
        end
        
        % m=n=0 时特殊处理
        if m == 0 && n == 0
            I_matrix(mi, ni) = 0;
            continue;
        end
        
        % 计算被积函数
        theta = m * dk1 + n * dk2;
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
        
        % 二重积分
        I_matrix(mi, ni) = trapz(k2, trapz(k1, integrand, 1));
    end
end

elapsed_time = toc;
fprintf('\n计算完成！总耗时: %.2f 秒\n', elapsed_time);

%% 2. 输出非零值
fprintf('\n========== 积分值 I(m,n) (m,n = -20 到 20) ==========\n');
fprintf('阈值: 1e-4\n\n');

% 找出非零值
[m_idx, n_idx] = find(abs(I_matrix) > 1e-4);
non_zero_count = length(m_idx);

fprintf('非零值总数: %d\n', non_zero_count);
fprintf('稀疏度: %.2f%%\n\n', non_zero_count/(length(m_vals)*length(n_vals))*100);

if non_zero_count > 0
    fprintf('m\tn\tI(m,n)\n');
    fprintf('----------------------------------------\n');
    for i = 1:non_zero_count
        m_val = m_vals(m_idx(i));
        n_val = n_vals(n_idx(i));
        I_val = I_matrix(m_idx(i), n_idx(i));
        fprintf('%3d\t%3d\t%.6f\n', m_val, n_val, I_val);
    end
end

%% 3. 按 |m|,|n| 分组输出
fprintf('\n========== 按 |m|, |n| 分组 ==========\n');
fprintf('|m|\t|n|\tI(m,n)\n');
fprintf('----------------------------------------\n');

% 收集分组数据
group_data = [];
for i = 1:non_zero_count
    m_abs = abs(m_vals(m_idx(i)));
    n_abs = abs(n_vals(n_idx(i)));
    I_val = I_matrix(m_idx(i), n_idx(i));
    group_data = [group_data; m_abs, n_abs, I_val];
end

% 按 |m| 排序输出
[~, sort_idx] = sort(group_data(:,1));
group_data = group_data(sort_idx, :);

for i = 1:size(group_data, 1)
    fprintf('%d\t%d\t%.6f\n', group_data(i,1), group_data(i,2), group_data(i,3));
end

%% 4. 绘制结果
figure('Position', [50, 50, 1400, 800]);

% 子图1: 热力图
subplot(2, 3, 1);
display_mat = I_matrix;
display_mat(abs(display_mat) < 1e-4) = 0;
imagesc(n_vals, m_vals, display_mat);
colorbar;
xlabel('n', 'FontSize', 12);
ylabel('m', 'FontSize', 12);
title('I(m,n) 热力图 (-20 到 20)', 'FontSize', 13);
set(gca, 'YDir', 'normal');
axis square;

% 子图2: 3D 山脉图
subplot(2, 3, 2);
[M, N] = meshgrid(m_vals, n_vals);
display_3d = I_matrix;
display_3d(abs(display_3d) < 1e-4) = NaN;
surf(M, N, display_3d', 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xlabel('m', 'FontSize', 11);
ylabel('n', 'FontSize', 11);
zlabel('I(m,n)', 'FontSize', 11);
title('3D 山脉结构', 'FontSize', 13);
colormap('jet');
colorbar;
grid on;
view(45, 30);
zlim([-1.2, 1.2]);

% 子图3: 沿 m 方向 (n=0)
subplot(2, 3, 3);
idx_n0 = find(n_vals == 0);
plot(m_vals, I_matrix(idx_n0, :), 'bo-', 'LineWidth', 2, 'MarkerSize', 4);
xlabel('m', 'FontSize', 12);
ylabel('I(m,0)', 'FontSize', 12);
title('沿 m 方向 (n=0)', 'FontSize', 13);
grid on;
xlim([-20, 20]);

% 子图4: 沿 n 方向 (m=0)
subplot(2, 3, 4);
idx_m0 = find(m_vals == 0);
plot(n_vals, I_matrix(idx_m0, :), 'ro-', 'LineWidth', 2, 'MarkerSize', 4);
xlabel('n', 'FontSize', 12);
ylabel('I(0,n)', 'FontSize', 12);
title('沿 n 方向 (m=0)', 'FontSize', 13);
grid on;
xlim([-20, 20]);

% 子图5: 次对角线 (|m|=|n|+1)
subplot(2, 3, 5);
diag_vals = zeros(length(m_vals), 1);
for mi = 1:length(m_vals)
    m = m_vals(mi);
    if m > 0
        n = m - 1;
        idx_n = find(n_vals == n);
        if ~isempty(idx_n)
            diag_vals(mi) = I_matrix(mi, idx_n);
        end
    end
end
plot(0:length(m_vals)-1, diag_vals, 'go-', 'LineWidth', 2, 'MarkerSize', 4);
xlabel('m', 'FontSize', 12);
ylabel('I(m,m-1)', 'FontSize', 12);
title('次对角线 (|m|=|n|+1)', 'FontSize', 13);
grid on;
xlim([0, 20]);

% 子图6: 统计信息
subplot(2, 3, 6);
text(0.1, 0.9, sprintf('范围: m,n = -20 到 20'), 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.75, sprintf('总点数: %d', length(m_vals)*length(n_vals)), 'FontSize', 11);
text(0.1, 0.6, sprintf('非零值: %d', non_zero_count), 'FontSize', 11);
text(0.1, 0.45, sprintf('稀疏度: %.2f%%', non_zero_count/(length(m_vals)*length(n_vals))*100), 'FontSize', 11);
if non_zero_count > 0
    text(0.1, 0.3, sprintf('最大值: %.6f', max(abs(I_matrix(:)))), 'FontSize', 11);
    text(0.1, 0.15, sprintf('最小值: %.6f', min(abs(I_matrix(abs(I_matrix)>1e-4)))), 'FontSize', 11);
end
axis off;
title('统计信息', 'FontSize', 13);

%% 5. 保存结果
save('Integral_large_range_20.mat', 'm_vals', 'n_vals', 'I_matrix');
fprintf('\n数据已保存到 Integral_large_range_20.mat\n');

%% 6. 输出到文本文件
fid = fopen('Integral_Values_20.txt', 'w');
fprintf(fid, '积分值 I(m,n) (m,n = -20 到 20)\n');
fprintf(fid, '========================================\n');
fprintf(fid, 'm\tn\tI(m,n)\n');
fprintf(fid, '----------------------------------------\n');

for i = 1:non_zero_count
    m_val = m_vals(m_idx(i));
    n_val = n_vals(n_idx(i));
    I_val = I_matrix(m_idx(i), n_idx(i));
    fprintf(fid, '%d\t%d\t%.6f\n', m_val, n_val, I_val);
end
fclose(fid);

fprintf('积分值已保存到 Integral_Values_20.txt\n');