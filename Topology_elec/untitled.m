% 计算二重积分 I(m,n) 并绘制沿m方向的子图（含数值拟合渐近背景）
% I(m,n) = ∫∫_Ω [1 - cos(mω₁ + nω₂)] / [2 - (cosω₁ + cosω₂)] dω₁ dω₂
% 渐近行为（m→∞）: I(m,n) ~ a(n)·ln(m) + b(n)

clear; clc; close all;

%% 1. 参数设置
m_vals = 0:0.5:80;        % m 取值范围
n_fixed = [0, 1, 2, 3, 5, 8, 10];   % 固定的 n 值
N_quad = 300;             % 积分网格精度

% 积分网格
omega1 = linspace(-pi, pi, N_quad);
omega2 = linspace(-pi, pi, N_quad);
[dW1, dW2] = meshgrid(omega1, omega2);

% 被积函数分母
denom = 2 - (cos(dW1) + cos(dW2));
denom(abs(denom) < 1e-12) = 1e-12;

%% 2. 数值计算 I(m,n)
colors = lines(length(n_fixed));
I_results = cell(length(n_fixed), 1);

fprintf('开始数值计算...\n');
for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    I_vals = zeros(size(m_vals));
    
    for mi = 1:length(m_vals)
        m = m_vals(mi);
        if m == 0
            numerator = 1 - cos(n * dW2);
        else
            numerator = 1 - cos(m * dW1 + n * dW2);
        end
        integrand = numerator ./ denom;
        I_vals(mi) = trapz(omega2, trapz(omega1, integrand, 1));
    end
    
    I_results{ni} = I_vals;
    fprintf('n = %d 完成\n', n);
end

%% 3. 拟合渐近系数（用大m区域）
% 拟合形式: I(m,n) = a(n) * ln(m) + b(n)
fit_a = zeros(size(n_fixed));
fit_b = zeros(size(n_fixed));
fit_r2 = zeros(size(n_fixed));

fprintf('\n拟合渐近系数（m > 50）:\n');
fprintf('n\ta\t\tb\t\tR²\n');

for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    
    % 取大m区域的数据
    idx_fit = m_vals > 50;
    m_fit = m_vals(idx_fit);
    I_fit = I_results{ni}(idx_fit);
    
    if length(m_fit) > 3
        % 线性拟合: I = a * ln(m) + b
        X = [log(m_fit)', ones(length(m_fit), 1)];
        coef = X \ I_fit';
        fit_a(ni) = coef(1);
        fit_b(ni) = coef(2);
        
        % 计算R²
        I_pred = X * coef;
        SS_res = sum((I_fit' - I_pred).^2);
        SS_tot = sum((I_fit' - mean(I_fit')).^2);
        fit_r2(ni) = 1 - SS_res / SS_tot;
        
        fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\n', n, fit_a(ni), fit_b(ni), fit_r2(ni));
    end
end

%% 4. 绘制主图：数值曲线 + 拟合渐近背景
figure('Name', '沿m方向子图（含渐近背景）', 'Position', [100, 100, 1100, 700]);
hold on;

% 4.1 绘制拟合渐近线作为背景（灰色虚线，只在大m区域）
m_asym_start = 30;  % 从m=30开始显示渐近线
m_asym = linspace(m_asym_start, max(m_vals)*1.1, 200);

for ni = 1:length(n_fixed)
    if fit_r2(ni) > 0.99  % 只显示拟合好的
        I_asym = fit_a(ni) * log(m_asym) + fit_b(ni);
        plot(m_asym, I_asym, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end
end

% 4.2 绘制数值曲线（彩色实线+标记）
for ni = 1:length(n_fixed)
    plot(m_vals, I_results{ni}, 'o-', 'LineWidth', 2, 'Color', colors(ni,:), ...
        'MarkerSize', 4, 'DisplayName', sprintf('n = %d (数值)', n_fixed(ni)));
end

xlabel('m', 'FontSize', 12);
ylabel('I(m,n)', 'FontSize', 12);
title('沿 m 方向变化（固定 n），灰色虚线为拟合渐近 I ~ a(n)·ln(m) + b(n)', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
xlim([0, max(m_vals)*1.05]);

% 添加文字说明
text(0.02, 0.98, '渐近形式: I(m,n) ~ a(n)·ln(m) + b(n)', ...
    'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', ...
    'VerticalAlignment', 'top');

% 添加竖线标记渐近区域的起始位置
ylim_current = ylim;
line([m_asym_start, m_asym_start], ylim_current, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
text(m_asym_start, ylim_current(2)*0.95, '渐近区域', ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'BackgroundColor', 'w');

hold off;

%% 5. 单独子图：n=0 时的详细对比
figure('Name', 'n=0 详细对比', 'Position', [150, 150, 900, 600]);

% 找到 n=0 的数据
idx_n0 = find(n_fixed == 0);
if ~isempty(idx_n0)
    I_n0 = I_results{idx_n0};
    
    % 数值曲线
    plot(m_vals, I_n0, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, ...
        'DisplayName', '数值积分');
    hold on;
    
    % 拟合渐近线（全范围）
    m_all = linspace(1, max(m_vals)*1.1, 200);
    I_asym_all = fit_a(idx_n0) * log(m_all) + fit_b(idx_n0);
    plot(m_all, I_asym_all, 'r--', 'LineWidth', 2, ...
        'DisplayName', sprintf('渐近: %.4f·ln(m) + %.4f', fit_a(idx_n0), fit_b(idx_n0)));
    
    % 标记拟合区域
    idx_fit = m_vals > 50;
    plot(m_vals(idx_fit), I_n0(idx_fit), 'go', 'MarkerSize', 8, ...
        'DisplayName', '拟合数据点 (m>50)');
    
    xlabel('m', 'FontSize', 12);
    ylabel('I(m,0)', 'FontSize', 12);
    title(sprintf('n=0: 数值积分 vs 渐近拟合 (R² = %.4f)', fit_r2(idx_n0)), 'FontSize', 13);
    legend('Location', 'northwest', 'FontSize', 10);
    grid on;
    xlim([0, max(m_vals)*1.05]);
    hold off;
end

%% 6. 绘制斜率 a(n) 和截距 b(n) 随 n 的变化
figure('Name', '渐近系数随n变化', 'Position', [200, 200, 1100, 400]);

subplot(1, 2, 1);
plot(n_fixed, fit_a, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('n', 'FontSize', 12);
ylabel('a(n)', 'FontSize', 12);
title('斜率 a(n) 随 n 变化', 'FontSize', 13);
grid on;

subplot(1, 2, 2);
plot(n_fixed, fit_b, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('n', 'FontSize', 12);
ylabel('b(n)', 'FontSize', 12);
title('截距 b(n) 随 n 变化', 'FontSize', 13);
grid on;

%% 7. 验证：I(m,n) - b(n) 与 ln(m) 的关系
figure('Name', '验证ln(m)关系', 'Position', [250, 250, 1100, 700]);
hold on;

for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    idx_positive = m_vals > 1;
    m_plot = m_vals(idx_positive);
    I_plot = I_results{ni}(idx_positive);
    
    % I - b 应该与 ln(m) 成正比
    I_minus_b = I_plot - fit_b(ni);
    
    plot(m_plot, I_minus_b, 'o-', 'LineWidth', 1.5, 'Color', colors(ni,:), ...
        'DisplayName', sprintf('n = %d', n));
end

xlabel('m', 'FontSize', 12);
ylabel('I(m,n) - b(n)', 'FontSize', 12);
title('验证: I - b(n) 应与 ln(m) 成正比', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'XScale', 'log');

% 添加参考线 (ln(m) 的斜率)
m_ref = [10, max(m_vals)];
ln_ref = log(m_ref);
for ni = 1:length(n_fixed)
    plot(m_ref, fit_a(ni) * ln_ref, '--', 'Color', colors(ni,:), ...
        'LineWidth', 1, 'HandleVisibility', 'off');
end

hold off;

%% 8. 输出最终结果
fprintf('\n最终渐近系数:\n');
fprintf('n\ta(n)\t\tb(n)\t\tR²\n');
for ni = 1:length(n_fixed)
    fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\n', ...
        n_fixed(ni), fit_a(ni), fit_b(ni), fit_r2(ni));
end

fprintf('\n结论: I(m,n) ~ a(n)·ln(m) + b(n), m→∞\n');
fprintf('其中 a(n) ≈ %.4f - %.4f·n (近似线性)\n', ...
    fit_a(1), (fit_a(1)-fit_a(end))/n_fixed(end));