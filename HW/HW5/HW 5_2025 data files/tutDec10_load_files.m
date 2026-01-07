clearvars, clc

M = readmatrix("mass_matrixMATLAB.xlsx");
K = readmatrix("stiffness_matrixMATLAB.xlsx");
C = readmatrix("damping_matrixMATLAB.xlsx");

format long
disp('Upper left sub-half-matrix of K')
K(1:4,1:4)
disp('Lower right sub-half-matrix of K')
K(5:8,5:8)
disp('Upper right sub-half-matrix of K')
K(1:4,5:8)

shouldbezero8times8 = abs(K - transpose(K));
disp('shouldbezero8times8:')
for ii = 1:8
    fprintf(' %12.10f',shouldbezero8times8(ii,:))
    fprintf('\n')
end
fprintf('\n')

disp('Upper left sub-half-matrix of M')
M(1:4,1:4)
disp('Lower right sub-half-matrix of M')
M(5:8,5:8)
disp('Upper right sub-half-matrix of M')
M(1:4,5:8)

shouldbezero8times8 = abs(M - transpose(M));
disp('shouldbezero8times8:')
for ii = 1:8
    fprintf(' %12.10f',shouldbezero8times8(ii,:))
    fprintf('\n')
end
fprintf('\n')

disp('Upper left sub-half-matrix of C')
C(1:4,1:4)
disp('Lower right sub-half-matrix of C')
C(5:8,5:8)
disp('Upper right sub-half-matrix of C')
C(1:4,5:8)

shouldbezero8times8 = abs(C - transpose(C));
disp('shouldbezero8times8:')
for ii = 1:8
    fprintf(' %12.10f',shouldbezero8times8(ii,:))
    fprintf('\n')
end
fprintf('\n')

% ========================================================
%  以下是添加的计算部分 (Calculation Part)
% ========================================================

% 1. 求解广义特征值问题 (Generalized Eigenvalue Problem)
% 求解方程: (K - lambda * M) * v = 0
% [V, D] = eig(K, M)
% V 是特征向量矩阵 (每一列是一个 Mode Shape)
% D 是对角矩阵，对角线上的元素是特征值 (omega^2)
[V, D] = eig(K, M);

% 2. 提取自然频率 (Natural Frequencies)
% 从对角矩阵提取特征值
eigenvalues = diag(D); 

% 特征值是频率的平方，所以要开根号得到 omega (rad/s)
omega_n = sqrt(eigenvalues);

% 3. 排序 (Sorting)
% 虽然 eig 通常会排序，但为了保险起见，我们手动按从小到大排序
[omega_sorted, sort_idx] = sort(omega_n);
V_sorted = V(:, sort_idx); % 同时要把对应的特征向量(Mode Shapes)也排好序

% 打印结果：自然频率
disp('Undamped Natural Frequencies (rad/s):');
disp(omega_sorted);

% 4. 计算作业要求的比值 (Ratio Calculation)
% 找到最低频率对应的 Mode Shape (排序后的第1列)
lowest_mode_shape = V_sorted(:, 1);

% 获取第1个位移 (Index 1) 和最后一个位移 (Index 8 或 end)
disp_first = lowest_mode_shape(1);
disp_last = lowest_mode_shape(end);

% 计算比值
ratio = disp_first / disp_last;

fprintf('\n=== Results for Homework 5 ===\n');
fprintf('Lowest Natural Frequency: %.5f rad/s\n', omega_sorted(1));
fprintf('Highest Natural Frequency: %.5f rad/s\n', omega_sorted(8));

fprintf('Displacement Ratio (1st / 8th): %.6f\n', ratio);