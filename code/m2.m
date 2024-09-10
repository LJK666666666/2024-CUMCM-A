clear
clc
close all

% 定义常数
delta_pho = 0.55; % 螺距，代表龙形路径的距离
a = delta_pho / (2 * pi); % 与螺距相关的常数
theta_1_0 = 16 * 2 * pi; % 龙头初始的角度位置（16圈）
v = 1; % 龙头的速度 (米/秒)
output = [1,2,52,102,152,202,224]; % 指定要输出结果的板凳索引

% 逐步增加时间，直到发生碰撞
% 依次细化时间步长，精确定位碰撞时刻

% for t=0:1000
% for t=300:0.1:1000
% for t=410:0.01:1000
% for t=412.47:0.001:1000
% for t=412.473:0.00001:1000
for t=412.47382:0.0000001:1000
    % 计算当前时间的龙头角度位置
    theta_1_t = cal_theta_t2(theta_1_0, v, t, a, 1);
    
    % 计算当前时间所有板凳的角度位置
    theta_list = cal_theta_all(theta_1_t, a);
    
    % 根据角度位置计算板凳的坐标
    p_list = cal_xy(theta_list, a);
    x_list = p_list(:, 1); % 获取所有板凳的 x 坐标
    y_list = p_list(:, 2); % 获取所有板凳的 y 坐标

    % 判断是否发生碰撞
    if judge_crash(theta_list, a) == 1
        break; % 如果发生碰撞，则退出循环
    end
end

disp('碰撞发生时刻: ')
disp(vpa(t));
disp(vpa(t,9));



% 使用已知的碰撞时间点进行进一步计算和输出
t = 412.473838; % 碰撞时刻的精确值
theta_1_t = cal_theta_t2(theta_1_0, v, t, a, 1); % 计算龙头在碰撞时刻的角度位置
theta_list = cal_theta_all(theta_1_t, a); % 计算所有板凳在碰撞时刻的角度位置

% 计算碰撞时刻所有板凳的坐标
p_list = cal_xy(theta_list, a);
x_list = p_list(:, 1); % 获取 x 坐标
y_list = p_list(:, 2); % 获取 y 坐标
v_list = cal_v_all(theta_list, v); % 计算各板凳的速度

result1 = [p_list,v_list];
writematrix(round(result1,6),'result2.xlsx','Range','B2');

result2 = p_list(output,:)';
disp('问题二位置结果: ');
disp(vpa(round(result2,6),8));

result3 = v_list(output)';
disp('问题二速度结果: ');
disp(vpa(round(result3,6),8));



% 绘制龙形的图像
figure;
scatter(x_list,y_list,10,'filled','MarkerFaceColor','blue');
hold on;
for i=1:223
    P = generate_dot(theta_list(i),theta_list(i+1),a);
    plot(P([1,2],1),P([1,2],2));
    hold on;
    plot(P([1,3],1),P([1,3],2));
    hold on;
    plot(P([2,4],1),P([2,4],2));
    hold on;
    plot(P([3,4],1),P([3,4],2));
    hold on;
end
axis equal;
