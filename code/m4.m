clear
clc
close all

global p_inter2;
output1 = [1,2,52,102,152,202,224];
output2 = [1,51,101,151,201];

% 定义螺距、半径、角度等参数
delta_pho = 1.7; % 盘入螺距
a = delta_pho / (2 * pi); % 与螺距相关的常数
r = 4.5; % 调头空间的半径
theta_inter1 = r / a; % 计算角度位置

% 计算三个关键点的坐标
p_inter1 = [a * theta_inter1 * cos(theta_inter1), a * theta_inter1 * sin(theta_inter1)]; % 第一个交点
p_inter2 = -p_inter1 / 3; % 第二个交点
p_inter3 = -p_inter1; % 第三个交点

% 计算切向方向和旋转方向
direction_tl = [cos(theta_inter1) - theta_inter1 * sin(theta_inter1), sin(theta_inter1) + theta_inter1 * cos(theta_inter1)];
direction_r = direction_tl * [0 1; -1 0];
theta_inter_ = atan(p_inter1(2) / p_inter1(1));
alpha = acos(abs(dot(direction_r, p_inter3) / norm(direction_r) / norm(p_inter3))); % 计算alpha角度
r1 = r / 3 * 2 / cos(alpha); % 计算第一个圆弧的半径
r2 = r1 / 2; % 计算第二个圆弧的半径

% 计算两个圆心的位置
p_center1 = p_inter1+r1.*[cos(theta_inter_-alpha),sin(theta_inter_-alpha)];
p_center2 = p_inter2*3/2-p_center1/2;


% 计算时间序列的角度位置
theta_t2_list_ = zeros(101,1);
theta_t2_list_(1) = theta_inter1;
v = 1;
for t=1:100
    theta_t2_list_(t+1) = cal_theta_t2(theta_inter1,v,t,a,-1);
end
theta_mat_ = zeros(224,101);
for i=1:101
    theta_mat_(:,i) = cal_theta_all(theta_t2_list_(i),a);
end

p_mat_ = zeros(224,101,2);
p_mat_(:,:,1) = a.*theta_mat_.*cos(theta_mat_);
p_mat_(:,:,2) = a.*theta_mat_.*sin(theta_mat_);
p1 = theta_mat_(output1,[1,51,101]);
p1_ = zeros(length(output1)*2,3);
for j=1:3
    p = cal_xy(p1(:,j),a);
    p = p';
    p1_(:,j) = p(:);
end

v_mat_ = cal_v_all(theta_mat_,v);

% p1_ : 0,-50,-100位置
% v_mat_(output,[1,51,101]) : 0,-50,-100速度


% 标记不同段的路径，1代表螺线部分，2代表第一个圆弧部分，3代表第二个圆弧部分，4代表逆向螺线部分
flag_mat = ones(224,101);
l1 = r1*(pi-alpha*2);  % 第一个圆弧的弧长
l2 = l1/2;  % 第二个圆弧的弧长
flag_mat(1,2:floor(l1/v+1)) = 2;
flag_mat(1,ceil(l1/v+1):floor((l1+l2)/v+1)) = 3;
flag_mat(1,ceil((l1+l2)/v+1):end) = 4;
p_mat = zeros(224,101,2);
theta_mat = zeros(224,101);
theta_mat(:,1) = theta_mat_(:,1);
p_mat(:,1,1) = a.*theta_mat(:,1).*cos(theta_mat(:,1));
p_mat(:,1,2) = a.*theta_mat(:,1).*sin(theta_mat(:,1));
% 处理第一个圆弧
for j = 2:floor(l1 / v + 1)
    lambda = theta_inter_ - alpha + pi - (j - 1) * v / r1;
    theta_mat(1, j) = lambda;
    p_mat(1, j, :) = p_center1 + r1 * [cos(lambda), sin(lambda)];
end

% 处理第二个圆弧
for j = ceil(l1 / v + 1):floor((l1 + l2) / v + 1)
    lambda = theta_inter_ + alpha - pi + ((j - 1) * v - l1) / r2;
    theta_mat(1, j) = lambda;
    p_mat(1, j, :) = p_center2 + r2 * [cos(lambda), sin(lambda)];
end

% 处理逆向螺线
for j = ceil((l1 + l2) / v + 1):size(flag_mat, 2)
    theta_mat(1, j) = cal_theta_t2(theta_inter1, v, ((j - 1) * v - l1 - l2) / v, a, -1);
    p_mat(1, j, :) = -a .* theta_mat(1, j) .* [cos(theta_mat(1, j)), sin(theta_mat(1, j))];
end
for i=2:224
    for j=2:101
        if flag_mat(i-1,j)==1
            flag_mat(i,j) = 1;
            theta_n1 = cal_theta_n1(theta_mat(i-1,j),cal_d(i-1),a);
            theta_mat(i,j) = theta_n1;
            p_mat(i,j,:) = a.*theta_n1.*[cos(theta_n1),sin(theta_n1)];
        end
        if flag_mat(i-1,j)==2
            if norm(reshape(p_mat(i-1,j,:),1,2,1)-p_inter1)>cal_d(i-1)
                flag_mat(i,j) = 2;
                lambda = theta_mat(i-1,j)+asin(cal_d(i-1)/2/r1)*2;
                theta_mat(i,j) = lambda;
                p_mat(i,j,:) = p_center1+r1*[cos(lambda),sin(lambda)];

            else
                flag_mat(i,j) = 1;
                pace = (cal_d(i-1)-norm(reshape(p_mat(i-1,j,:),1,2,1)-p_inter1))/r/2;
                if pace<1e-9
                    pace=1e-9;
                end
                theta_estimate = theta_inter1;
                while pace>1e-10
                    theta_estimate = theta_estimate+pace;
                    if norm(reshape(p_mat(i-1,j,:),1,2,1)-a*theta_estimate*[cos(theta_estimate),sin(theta_estimate)])>cal_d(i-1)
                        theta_estimate = theta_estimate-pace;
                        pace = pace/2;
                    end
                end
                theta_mat(i,j) = theta_estimate;
                p_mat(i,j,:) = a*theta_estimate*[cos(theta_estimate),sin(theta_estimate)];

            end
        end
        if flag_mat(i-1,j)==3
            if norm(reshape(p_mat(i-1,j,:),1,2,1)-p_inter2)>cal_d(i-1)
                flag_mat(i,j) = 3;
                lambda = theta_mat(i-1,j)-asin(cal_d(i-1)/2/r2)*2;
                theta_mat(i,j) = lambda;
                p_mat(i,j,:) = p_center2+r2*[cos(lambda),sin(lambda)];

            else
                flag_mat(i,j) = 2;
                pace = (cal_d(i-1)-norm(reshape(p_mat(i-1,j,:),1,2,1)-p_inter2))/r1/2;
                if pace<1e-9
                    pace=1e-9;
                end
                theta_estimate = theta_inter_+alpha;
                while pace>1e-10
                    theta_estimate = theta_estimate+pace;
                    if norm(reshape(p_mat(i-1,j,:),1,2,1)-(p_center1+r1*[cos(theta_estimate),sin(theta_estimate)]))>cal_d(i-1)
                        theta_estimate = theta_estimate-pace;
                        pace = pace/2;
                    end
                end
                theta_mat(i,j) = theta_estimate;
                p_mat(i,j,:) = p_center1+r1*[cos(theta_estimate),sin(theta_estimate)];

            end

        end
        if flag_mat(i-1,j)==4
            if flag_mat(i,j-1)==4 | norm(reshape(p_mat(i-1,j,:),1,2,1)-p_inter3)>cal_d(i-1)
                flag_mat(i,j) = 4;
                theta_n1 = cal_theta_n1_reverse(theta_mat(i-1,j),cal_d(i-1),a);
                theta_mat(i,j) = theta_n1;
                p_mat(i,j,:) = -a.*theta_n1.*[cos(theta_n1),sin(theta_n1)];

            else
                flag_mat(i,j) = 3;
                pace = (cal_d(i-1)-norm(reshape(p_mat(i-1,j,:),1,2,1)-p_inter3))/r2/2;
                if pace<1e-9
                    pace=1e-9;
                end
                theta_estimate = theta_inter_-alpha;
                while pace>1e-10
                    theta_estimate = theta_estimate-pace;
                    if norm(reshape(p_mat(i-1,j,:),1,2,1)-(p_center2+r2*[cos(theta_estimate),sin(theta_estimate)]))>cal_d(i-1)
                        theta_estimate = theta_estimate+pace;
                        pace = pace/2;
                    end
                end
                theta_mat(i,j) = theta_estimate;
                p_mat(i,j,:) = p_center2+r2*[cos(theta_estimate),sin(theta_estimate)];

            end
        end
    end     
end
result01 = p_mat(output1,:,:);
result02 = zeros(length(output1)*2,101);
for i=1:101
     result03 = result01(:,i,:);
     result03 = permute(result03,[3 1 2]);
     result03 = result03(:);
     result02(:,i) = result03;
end


n0_mat = p_mat(2:224,:,:)-p_mat(1:223,:,:);
n1_mat = cal_n1_all(theta_mat,flag_mat);
v_mat = zeros(224,101);
v_mat(1,:) = v;
for i=1:101
    for j=2:224
        v_mat(j,i) = cal_v2(v_mat(j-1,i),n0_mat(j-1,i,:),n1_mat(j-1,i,:),n1_mat(j,i,:));
    end
end

v_mat_all = [v_mat_(:,101:-1:2),v_mat];
p_mat_all = [p_mat_(:,101:-1:2,:),p_mat];

result1 = permute(p_mat_all,[3,1,2]);
result1 = reshape(result1,448,201,1);
writematrix(round(result1,6),'result4.xlsx','Sheet','位置','Range','B2');

result2 = v_mat_all;
writematrix(round(result2,6),'result4.xlsx','Sheet','速度','Range','B2');


% 绘制速度变化图
figure;
mesh(-100:100,1:224,v_mat_all);
xlabel('时刻/t');
ylabel('把手序号');
zlabel('速度/(m·s^{-1})');

% 绘制龙形的末尾状态图
figure;
scatter(p_mat_all(:,end,1),p_mat_all(:,end,2),10,'filled','MarkerFaceColor','blue');
hold on;
for i=1:223
    P = generate_dot_(p_mat_all(i,end,:),p_mat_all(i+1,end,:));
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



% 绘制调头区域与龙形的路径
figure;

width1 = 1.5;
width2 = 2;

R_ = r;                    % 半径
center_ = [0,0];          % 圆心坐标
theta_start_ = 0;          % 起始角度（弧度）
theta_end_ = pi*2;         % 结束角度（弧度）
theta_ = linspace(theta_start_, theta_end_, 100);  % 角度范围
x_ = center_(1) + R_ * cos(theta_);  % 参数化方程：x(t) = x0 + R * cos(t)
y_ = center_(2) + R_ * sin(theta_);  % 参数化方程：y(t) = y0 + R * sin(t)
plot(x_, y_, 'LineWidth', width2);  % 绘制圆弧
hold on;

R_ = r1;                    % 半径
center_ = p_center1;          % 圆心坐标
theta_start_ = theta_inter_+alpha;          % 起始角度（弧度）
theta_end_ = pi+theta_inter_-alpha;         % 结束角度（弧度）
theta_ = linspace(theta_start_, theta_end_, 100);  % 角度范围
x_ = center_(1) + R_ * cos(theta_);  % 参数化方程：x(t) = x0 + R * cos(t)
y_ = center_(2) + R_ * sin(theta_);  % 参数化方程：y(t) = y0 + R * sin(t)
plot(x_, y_, 'black', 'LineWidth', width1);  % 绘制圆弧
hold on;

R_ = r2;                    % 半径
center_ = p_center2;          % 圆心坐标
theta_start_ = theta_inter_+alpha-pi;          % 起始角度（弧度）
theta_end_ = theta_inter_-alpha;         % 结束角度（弧度）
theta_ = linspace(theta_start_, theta_end_, 100);  % 角度范围
x_ = center_(1) + R_ * cos(theta_);  % 参数化方程：x(t) = x0 + R * cos(t)
y_ = center_(2) + R_ * sin(theta_);  % 参数化方程：y(t) = y0 + R * sin(t)
plot(x_, y_, 'black', 'LineWidth', width1);  % 绘制圆弧
hold on;

plot([p_center1(1),p_inter1(1)], [p_center1(2),p_inter1(2)], '--black', 'LineWidth', width1);
hold on;

plot([p_center1(1),p_inter2(1)], [p_center1(2),p_inter2(2)], '--black', 'LineWidth', width1);
hold on;

plot([p_center2(1),p_inter2(1)], [p_center2(2),p_inter2(2)], '--black', 'LineWidth', width1);
hold on;

plot([p_center2(1),p_inter3(1)], [p_center2(2),p_inter3(2)], '--black', 'LineWidth', width1);
hold on;

theta_ = linspace(theta_inter1, 7*pi, 1000);
x_ = a.*theta_.*cos(theta_);
y_ = a.*theta_.*sin(theta_);
plot(x_, y_, 'black', 'LineWidth', width1);  % 绘制圆弧
hold on;
plot(-x_, -y_, 'black', 'LineWidth', width1);  % 绘制圆弧
hold on;

dot_mat_ = [p_center1;p_center2;p_inter1;p_inter2;p_inter3];
scatter(dot_mat_(:,1),dot_mat_(:,2),15,'filled','black');

axis equal;



figure;
[X, Y] = meshgrid(-100:100,1:224);
scatter3(X(:), Y(:), v_mat_all(:), 10, v_mat_all(:), 'filled');

result3 = [p1_(:,[3 2]),result02(:,[1 51 101])];
disp('问题4位置结果: ')
disp(vpa(round(result3,6),8));

result4 = v_mat_all(output1,output2);
disp('问题4速度结果: ');
disp(vpa(round(result4,6),8));

scatter_animation(p_mat_all,r);
animation(p_mat_all,r);