%% 计算 V_A^n 和 V_B^n 的脚本
% 根据特征方程 lambda^2 + (C1/C2)*lambda + 1 = 0
% 递推关系：V_A^n = -(C1/C2)*V_A^(n-1) + V_B^(n-1)
%           V_B^n = -V_A^(n-1)
% 边界条件：V_A^0 给定，V_A^N = 0
% 适用于 C1 > 2*C2 的情况（实数特征根，振荡行为）

clear; clc; close all;

%% 参数设置
C1 = 3;      % 必须满足 C1 > 2*C2
C2 = 1;      % C2 > 0
N = 10;      % 可以取较大的值
V_A0 = 1;    % 初始值 V_A^0

%% 检查参数条件
if C1 <= 2*C2
    error('错误：需要 C1 > 2*C2 才能保证实数特征根，当前 C1 = %.2f, 2*C2 = %.2f', C1, 2*C2);
end

%% 计算特征根
lambda_plus = -C1/(2*C2) + sqrt((C1/(2*C2))^2 - 1);
lambda_minus = -C1/(2*C2) - sqrt((C1/(2*C2))^2 - 1);

fprintf('特征根：lambda_+ = %.6f, lambda_- = %.6f\n', lambda_plus, lambda_minus);

%% 计算系数 alpha 和 beta（公式42）
alpha = V_A0 / (1 - (lambda_plus/lambda_minus)^N);
beta = -V_A0 / (1 - (lambda_minus/lambda_plus)^N);

fprintf('系数：alpha = %.6e, beta = %.6e\n', alpha, beta);

%% 计算 V_A^n 序列（公式43）
n = (0:N)';
V_A_n = alpha * lambda_plus.^n + beta * lambda_minus.^n;

% 验证边界条件 V_A^N = 0
fprintf('V_A^N 的值（应接近0）：%.6e\n', V_A_n(end));

%% 计算 V_B^n 序列
% 根据递推关系：
% V_B^n = -V_A^(n-1)，其中 n = 1, 2, ..., N-1
% 由 V_A^1 = -(C1/C2)*V_A^0 + V_B^0 反推 V_B^0

V_B_n = zeros(N+1, 1);  % 预分配，索引从0开始

% 首先计算 V_B^0（从第一个递推式反推）
V_B_n(1) = V_A_n(2) + (C1/C2)*V_A_n(1);  % V_B^0 = V_A^1 + (C1/C2)*V_A^0

% 然后根据 V_B^n = -V_A^(n-1) 计算 n = 1 到 N-1
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
    V_A_recur = -(C1/C2)*V_A_n(n_idx) + V_B_n(n_idx);  % 使用 V_B^(n-1)
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
title(sprintf('V_A^n (C_1 = %.2f, C_2 = %.2f, N = %d)', C1, C2, N));
grid on;

% 子图2：V_B^n
subplot(2,1,2);
plot(0:N, V_B_n, 'r-', 'LineWidth', 1.5);
hold on;
plot(0:N, V_B_n, 'ro', 'MarkerSize', 6);
xlabel('n');
ylabel('V_B^n');
title(sprintf('V_B^n (C_1 = %.2f, C_2 = %.2f, N = %d)', C1, C2, N));
grid on;

% 子图3：同时显示
figure;
plot(0:N, V_A_n, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(0:N, V_B_n, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('n');
ylabel('值');
title(sprintf('V_A^n 和 V_B^n 对比 (C_1 = %.2f, C_2 = %.2f, N = %d)', C1, C2, N));
legend('V_A^n', 'V_B^n', 'Location', 'best');
grid on;
hold off;

%% 保存数据到表格
T = table((0:N)', V_A_n, V_B_n, 'VariableNames', {'n', 'V_A^n', 'V_B^n'});
disp(T);
% writetable(T, 'V_A_V_B_results.csv');  % 取消注释可保存为CSV文件

fprintf('\n计算完成！\n');