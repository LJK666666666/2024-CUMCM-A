clear
clc
close all

% 逐步减小螺距，寻找最小螺距
% 逐步细化搜索步长，提高搜索精度

% init = 0.55; pace = -0.01;
% init = 0.46; pace = -0.001;
% init = 0.451; pace = -0.0001;
% init = 0.4504; pace = -0.00001;
% init = 0.45034; pace = -0.000001;
init = 0.450338; pace = -0.0000001;
for delta_pho = init:pace:0.1
    a = delta_pho / (2 * pi); % 计算与螺距相关的常数
    theta_1_0 = 16 * 2 * pi; % 龙头初始角度位置（16圈）
    v = 1; % 龙头速度 (米/秒)
    output = [1, 2, 52, 102, 152, 202, 224]; % 指定要输出结果的板凳索引
    t = 0; % 时间初始值
    flag = 0; % 碰撞标志，初始为0（无碰撞）
    
    % 循环计算直到龙头达到或接近调头区域边界
    while a * cal_theta_t2(theta_1_0, v, t, a, 1) >= 4.5
        theta_1_t = cal_theta_t2(theta_1_0, v, t, a, 1); % 计算当前时间龙头的角度位置
        theta_list = cal_theta_all(theta_1_t, a); % 计算当前时间所有板凳的角度位置
        
        p_list = cal_xy(theta_list, a); % 计算所有板凳的坐标
        x_list = p_list(:, 1); % 获取所有板凳的 x 坐标
        y_list = p_list(:, 2); % 获取所有板凳的 y 坐标
        
        if judge_crash(theta_list, a) == 1 % 判断是否发生碰撞
            flag = 1; % 发生碰撞，标志设为1
            break; % 退出循环
        end
        
        t = t + 1; % 增加时间步长
        if a * cal_theta_t2(theta_1_0, v, t, a, 1) < 5
            t = t - 0.8; % 调整时间步长以避免跳过边界
        end
        if a * cal_theta_t2(theta_1_0, v, t, a, 1) < 4.6
            t = t - 0.19; % 进一步调整时间步长
        end
    end
    
    if flag == 1 % 如果发生了碰撞，终止搜索
        break;
    end
end

disp('最小螺距: ');
disp(vpa(delta_pho-pace));
disp(vpa(ceil((delta_pho-pace)*1e6)*1e-6,6));



% 绘制最终的龙形位置
delta_pho = 0.450338;

a = delta_pho/2/pi;
theta_1_0 = 16*2*pi;
v = 1;
output = [1,2,52,102,152,202,224];
t = 0;
while a*cal_theta_t2(theta_1_0,v,t,a,1)>=4.5
    theta_1_t = cal_theta_t2(theta_1_0,v,t,a,1);
    theta_list = cal_theta_all(theta_1_t,a);
    
    p_list = cal_xy(theta_list,a);
    x_list = p_list(:,1);
    y_list = p_list(:,2);
    
    t = t+1;
    if a*cal_theta_t2(theta_1_0,v,t,a,1)<5
        t = t-0.8;  % 调整时间步长以避免跳过边界

    end
    if a*cal_theta_t2(theta_1_0,v,t,a,1)<4.6
        t = t-0.19;  % 进一步调整时间步长
    end
end



% 绘制龙形位置图
figure;
scatter(x_list,y_list,20,'filled','MarkerFaceColor','blue');
hold on;
viscircles([0 0],4.5);
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
