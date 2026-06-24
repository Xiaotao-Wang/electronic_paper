%% 计算 V_A^n 和 V_B^n 的脚本（包含两种情况）
% 情况1：C1 > 2*C2（过阻尼，实数特征根）
% 情况2：C1 = 2*C2（临界阻尼，重根 lambda = -1）
% 递推关系：V_A^n = -(C1/C2)*V_A^(n-1) + V_B^(n-1)
%           V_B^n = -V_A^(n-1)
% 边界条件：V_A^0 给定，V_A^N = 0

clear; clc; close all;

%% 参数设置
% 选择计算模式：'overdamped' 或 'critical'
mode = 'critical';  % 改为 'overdamped' 可计算过阻尼情况

C1 = 2;      % 对于临界阻尼，C1 = 2*C2
C2 = 1;      % C2 > 0
N = 30;      % 可以取较大的值
V_A0 = 1;    % 初始值 V_A^0

%% 根据模式计算 V_A^n
fprintf('========== 计算模式：%s ==========\n', mode);

if strcmp(mode, 'overdamped')
    % 情况(a)：C1 > 2*C2（过阻尼）
    if C1 <= 2*C2
        error('错误：过阻尼需要 C1 > 2*C2，当前 C1 = %.2f, 2*C2 = %.2f', C1, 2*C2);
    end
    
    % 计算特征根
    lambda_plus = -C1/(2*C2) + sqrt((C1/(2*C2))^2 - 1);
    lambda_minus = -C1/(2*C2) - sqrt((C1/(2*C2))^2 - 1);
    
    fprintf('特征根：lambda_+ = %.6f, lambda_- = %.6f\n', lambda_plus, lambda_minus);
    
    % 计算系数 alpha 和 beta
    alpha = V_A0 / (1 - (lambda_plus/lambda_minus)^N);
    beta = -V_A0 / (1 - (lambda_minus/lambda_plus)^N);
    
    fprintf('系数：alpha = %.6e, beta = %.6e\n', alpha, beta);
    
    % 计算 V_A^n
    n = (0:N)';
    V_A_n = alpha * lambda_plus.^n + beta * lambda_minus.^n;
    
elseif strcmp(mode, 'critical')
    % 情况(b)：C1 = 2*C2（临界阻尼）
    if abs(C1 - 2*C2) > 1e-10
        warning('临界阻尼需要 C1 = 2*C2，当前 C1 = %.2f, 2*C2 = %.2f', C1, 2*C2);
    end
    
    fprintf('临界阻尼：lambda = -1 (重根)\n');
    
    % 计算 V_A^n
    n = (0:N)';
    V_A_n = V_A0 * (1 - n/N) .* (-1).^n;
    
else
    error('未知模式，请选择 ''overdamped'' 或 ''critical''');
end

% 验证边界条件 V_A^N = 0
fprintf('V_A^N 的值（应接近0）：%.6e\n', V_A_n(end));

%% 计算 V_B^n 序列（根据递推关系）
% 递推关系：V_B^n = -V_A^(n-1)，其中 n = 1, 2, ..., N-1
% V_B^0 由 V_A^1 = -(C1/C2)*V_A^0 + V_B^0 反推

V_B_n = zeros(N+1, 1);  % 预分配，索引从0开始

% 计算 V_B^0（从第一个递推式反推）
V_B_n(1) = V_A_n(2) + (C1/C2)*V_A_n(1);  % V_B^0 = V_A^1 + (C1/C2)*V_A^0

% 根据 V_B^n = -V_A^(n-1) 计算 n = 1 到 N-1
for n_idx = 1:N-1
    V_B_n(n_idx+1) = -V_A_n(n_idx);  % V_B^n = -V_A^(n-1)
end

% 延伸计算 V_B^N（虽然递推式中没有定义，但可以计算）
V_B_n(N+1) = -V_A_n(N);  % V_B^N = -V_A^(N-1)

%% 验证递推关系 V_A^n = -(C1/C2)*V_A^(n-1) + V_B^(n-1)
fprintf('\n验证递推关系：\n');
fprintf('  n\t V_A^n(计算)\t V_A^n(递推)\t 误差\n');
fprintf('------------------------------------------------------\n');
for n_idx = 1:N
    V_A_recur = -(C1/C2)*V_A_n(n_idx) + V_B_n(n_idx);
    error_val = abs(V_A_n(n_idx+1) - V_A_recur);
    fprintf('  %d\t %.6f\t %.6f\t %.2e\n', n_idx, V_A_n(n_idx+1), V_A_recur, error_val);
end

%% 输出详细结果
fprintf('\n========== 计算结果 ==========\n');
fprintf('  n\t V_A^n\t\t V_B^n\n');
fprintf('----------------------------------------\n');
for i = 1:N+1
    if i == 1
        fprintf('  %d\t %.6f\t %.6f (反推)\n', i-1, V_A_n(i), V_B_n(i));
    elseif i == N+1
        fprintf('  %d\t %.6f\t %.6f (延伸)\n', i-1, V_A_n(i), V_B_n(i));
    else
        fprintf('  %d\t %.6f\t %.6f\n', i-1, V_A_n(i), V_B_n(i));
    end
end

%% 绘图
figure;

% 子图1：V_A^n
subplot(2,1,1);
plot(0:N, V_A_n, 'b-', 'LineWidth', 1.5);
hold on;
plot(0:N, V_A_n, 'bo', 'MarkerSize', 6);
xlabel('n');
ylabel('V_A^n');
title(sprintf('V_A^n (模式: %s, C_1 = %.2f, C_2 = %.2f, N = %d)', mode, C1, C2, N));
grid on;

% 子图2：V_B^n
subplot(2,1,2);
plot(0:N, V_B_n, 'r-', 'LineWidth', 1.5);
hold on;
plot(0:N, V_B_n, 'ro', 'MarkerSize', 6);
xlabel('n');
ylabel('V_B^n');
title(sprintf('V_B^n (模式: %s, C_1 = %.2f, C_2 = %.2f, N = %d)', mode, C1, C2, N));
grid on;

% 同时显示
figure;
plot(0:N, V_A_n, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(0:N, V_B_n, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('n');
ylabel('电压值');
title(sprintf('V_A^n 和 V_B^n 对比 (模式: %s, N = %d)', mode, N));
legend('V_A^n', 'V_B^n', 'Location', 'best');
grid on;
hold off;

%% 保存数据到表格
T = table((0:N)', V_A_n, V_B_n, 'VariableNames', {'n', 'V_A^n', 'V_B^n'});
disp(T);
% writetable(T, 'V_A_V_B_results.csv');  % 取消注释可保存为CSV文件

fprintf('\n计算完成！\n');