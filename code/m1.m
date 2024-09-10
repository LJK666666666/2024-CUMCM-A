clear
clc
close all

% 定义常数
delta_pho = 0.55; % 螺距，代表龙形路径的距离
a = delta_pho / (2 * pi); % 与螺距相关的常数
theta_1_0 = 16 * 2 * pi; % 龙头初始的角度位置（16圈）
v = 1; % 龙头的速度 (米/秒)
p = zeros(224, 301, 2); % 用于存储每节板凳把手的坐标的位置矩阵
theta_mat = zeros(224, 301); % 用于存储每节板凳的角度位置的矩阵
output1 = [1,2,52,102,152,202,224]; % 指定要输出结果的板凳索引
output2 = 1:60:301; % 指定要输出结果的时间点索引

% 计算龙头在每个时间点的角度位置
theta_t_list = zeros(1, 301); % 初始化角度位置列表
theta_t_list(1) = theta_1_0; % 设置初始角度位置
for i = 2:301
    theta_t_list(i) = cal_theta_t2(theta_1_0, v, i - 1, a, 1); % 计算每秒的角度位置
end

% 计算各个板凳的位置并存储在矩阵中
for i = 1:301
    theta_n_list = cal_theta_all(theta_t_list(i), a); % 计算当前时间下所有板凳的角度位置
    theta_mat(:, i) = theta_n_list; % 将角度位置存储在矩阵中
    p_n_list = cal_xy(theta_n_list, a); % 计算当前时间下所有板凳的坐标
    p(:, i, :) = p_n_list; % 将坐标存储在位置矩阵中
end

% 提取指定板凳在指定时间点的位置数据
result1 = p(output1, output2, :); % 获取指定板凳和时间点的位置
result1 = permute(result1, [3, 1, 2]); % 调整维度以便处理
result1 = reshape(result1, length(output1) * 2, length(output2), 1); % 重塑为二维矩阵
disp('问题一位置结果: ');
disp(vpa(round(result1, 6), 8)); % 显示位置数据，保留6位小数

% 计算速度并存储在矩阵中
v_mat = cal_v_all(theta_mat, v); % 计算所有板凳在所有时间点的速度
result2 = v_mat(output1, output2); % 提取指定板凳和时间点的速度
disp('问题一速度结果: ');
disp(vpa(round(result2, 6), 7)); % 显示速度数据，保留6位小数

% 将位置数据保存到Excel文件
result3 = permute(p, [3, 1, 2]); % 调整维度以便处理
result3 = reshape(result3, 448, 301, 1); % 重塑为二维矩阵
writematrix(round(result3, 6), 'result1.xlsx', 'Sheet', '位置', 'Range', 'B2'); % 保存到Excel文件

% 将速度数据保存到Excel文件
result4 = v_mat; % 使用速度矩阵直接保存
writematrix(round(result4, 6), 'result1.xlsx', 'Sheet', '速度', 'Range', 'B2'); % 保存到Excel文件

% 绘制龙形初始配置的图像
figure;
scatter(p(:, 1, 1), p(:, 1, 2), 10, 'filled', 'MarkerFaceColor', 'blue'); 
hold on;
for i = 1:223
    P = generate_dot(theta_mat(i, 1), theta_mat(i + 1, 1), a); 
    plot(P([1, 2], 1), P([1, 2], 2)); 
    hold on;
    plot(P([1, 3], 1), P([1, 3], 2)); 
    hold on;
    plot(P([2, 4], 1), P([2, 4], 2)); 
    hold on;
    plot(P([3, 4], 1), P([3, 4], 2)); 
    hold on;
end
axis equal; % 设置坐标轴比例相等



scatter_animation(p,0);
animation(p,0);
