% 计算二重积分 I(m,n) 并拟合渐近行为
clear; clc; close all;

%% 1. 参数设置
m_vals = 0:1:80;          % m 取值范围
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

%% 3. 拟合渐近行为
% 尝试不同形式的渐近函数，看哪个拟合最好
fprintf('\n渐近行为拟合（m > 50）:\n');
fprintf('n\tln(m)\t\tln(m)/√m\t\t1/√m\t\t常数\n');

for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    
    % 取大m区域的数据
    idx_fit = m_vals > 50;
    m_fit = m_vals(idx_fit);
    I_fit = I_results{ni}(idx_fit);
    
    if length(m_fit) > 5
        % 尝试不同的渐近形式
        % 1. I ~ a * ln(m) + b
        X1 = [log(m_fit)', ones(length(m_fit), 1)];
        coef1 = X1 \ I_fit';
        r2_1 = 1 - sum((I_fit' - X1*coef1).^2) / sum((I_fit' - mean(I_fit')).^2);
        
        % 2. I ~ a * ln(m)/√m + b
        X2 = [log(m_fit)'./sqrt(m_fit)', ones(length(m_fit), 1)];
        coef2 = X2 \ I_fit';
        r2_2 = 1 - sum((I_fit' - X2*coef2).^2) / sum((I_fit' - mean(I_fit')).^2);
        
        % 3. I ~ a / √m + b
        X3 = [1./sqrt(m_fit)', ones(length(m_fit), 1)];
        coef3 = X3 \ I_fit';
        r2_3 = 1 - sum((I_fit' - X3*coef3).^2) / sum((I_fit' - mean(I_fit')).^2);
        
        % 4. I ~ a * √m + b (检查是否真的递增)
        X4 = [sqrt(m_fit)', ones(length(m_fit), 1)];
        coef4 = X4 \ I_fit';
        r2_4 = 1 - sum((I_fit' - X4*coef4).^2) / sum((I_fit' - mean(I_fit')).^2);
        
        fprintf('%d\tR²=%.4f\t\tR²=%.4f\t\tR²=%.4f\t\tR²=%.4f\n', ...
            n, r2_1, r2_2, r2_3, r2_4);
        
        % 存储最佳拟合参数（假设ln(m)拟合最好）
        fit_coef{ni} = coef1;
        fit_type{ni} = 'ln(m)';
    end
end

%% 4. 绘制数值曲线 + 拟合渐近线
figure('Name', '数值结果与拟合渐近线', 'Position', [100, 100, 1100, 700]);
hold on;

for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    
    % 数值曲线
    plot(m_vals, I_results{ni}, 'o-', 'LineWidth', 2, 'Color', colors(ni,:), ...
        'MarkerSize', 4, 'DisplayName', sprintf('n = %d (数值)', n));
    
    % 拟合渐近线（使用ln(m)形式）
    idx_fit = m_vals > 50;
    m_fit = m_vals(idx_fit);
    I_fit = I_results{ni}(idx_fit);
    
    if length(m_fit) > 5
        X = [log(m_fit)', ones(length(m_fit), 1)];
        coef = X \ I_fit';
        
        % 绘制拟合线
        m_asym = linspace(50, max(m_vals)*1.1, 100);
        I_asym = coef(1) * log(m_asym) + coef(2);
        plot(m_asym, I_asym, '--', 'Color', colors(ni,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('n = %d (拟合: %.3f ln(m) %+.3f)', n, coef(1), coef(2)));
    end
end

xlabel('m', 'FontSize', 12);
ylabel('I(m,n)', 'FontSize', 12);
title('数值结果 vs 拟合渐近线（ln(m)形式）', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
xlim([0, max(m_vals)*1.05]);

hold off;

%% 5. 半对数坐标（检查ln(m)行为）
figure('Name', '半对数坐标验证', 'Position', [150, 150, 1100, 700]);
hold on;

for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    idx_positive = m_vals > 1;
    semilogx(m_vals(idx_positive), I_results{ni}(idx_positive), 'o-', ...
        'LineWidth', 2, 'Color', colors(ni,:), ...
        'MarkerSize', 4, 'DisplayName', sprintf('n = %d', n));
end

xlabel('m (对数坐标)', 'FontSize', 12);
ylabel('I(m,n)', 'FontSize', 12);
title('半对数坐标：如果渐近是ln(m)，应显示为直线', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10);
grid on;

hold off;

%% 6. 检查增长率
figure('Name', '增长率分析', 'Position', [200, 200, 1100, 700]);
hold on;

for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    
    % 计算差分增长率
    idx_positive = m_vals > 1;
    m_plot = m_vals(idx_positive);
    I_plot = I_results{ni}(idx_positive);
    
    % 计算 dI/dm
    dI = diff(I_plot) ./ diff(m_plot);
    m_mid = (m_plot(1:end-1) + m_plot(2:end)) / 2;
    
    plot(m_mid, dI, 'o-', 'LineWidth', 1.5, 'Color', colors(ni,:), ...
        'DisplayName', sprintf('n = %d', n));
end

xlabel('m', 'FontSize', 12);
ylabel('dI/dm (增长率)', 'FontSize', 12);
title('增长率分析：如果I ~ ln(m)，则dI/dm ~ 1/m', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');

hold off;

%% 7. 输出拟合结果
fprintf('\n最佳拟合参数（I ≈ a*ln(m) + b）:\n');
fprintf('n\ta\t\tb\t\tR²\n');
for ni = 1:length(n_fixed)
    n = n_fixed(ni);
    idx_fit = m_vals > 50;
    m_fit = m_vals(idx_fit);
    I_fit = I_results{ni}(idx_fit);
    
    if length(m_fit) > 5
        X = [log(m_fit)', ones(length(m_fit), 1)];
        coef = X \ I_fit';
        r2 = 1 - sum((I_fit' - X*coef).^2) / sum((I_fit' - mean(I_fit')).^2);
        fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\n', n, coef(1), coef(2), r2);
    end
end