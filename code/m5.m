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

v = 1;

time_inter = 0.001;
% time_inter = 0.5; 龙头最大行进速度: 1.247010m/s 
% time_inter = 0.1; 龙头最大行进速度: 1.247010m/s 
% time_inter = 0.01; 龙头最大行进速度: 1.246266m/s 
% time_inter = 0.001; 龙头最大行进速度: 1.246266m/s 

v = v*time_inter;

% 标记不同段的路径，1代表螺线部分，2代表第一个圆弧部分，3代表第二个圆弧部分，4代表逆向螺线部分
flag_mat = ones(224,100/time_inter+1);
l1 = r1*(pi-alpha*2);  % 第一个圆弧的弧长
l2 = l1/2;  % 第二个圆弧的弧长
flag_mat(1,2:floor(l1/v+1)) = 2;
flag_mat(1,ceil(l1/v+1):floor((l1+l2)/v+1)) = 3;
flag_mat(1,ceil((l1+l2)/v+1):end) = 4;
p_mat = zeros(224,100/time_inter+1,2);
theta_mat = zeros(224,100/time_inter+1);
theta_mat(:,1) = cal_theta_all(theta_inter1,a);
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
    for j=2:100/time_inter+1
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



n0_mat = p_mat(2:224,:,:)-p_mat(1:223,:,:);
n1_mat = cal_n1_all(theta_mat,flag_mat);
v_mat = zeros(224,100/time_inter+1);
v_mat(1,:) = v;
for i=1:100/time_inter+1
    for j=2:224
        v_mat(j,i) = cal_v2(v_mat(j-1,i),n0_mat(j-1,i,:),n1_mat(j-1,i,:),n1_mat(j,i,:));
    end
end

figure;
mesh(0:time_inter:100,1:224,v_mat);
xlabel('时刻/t');
ylabel('把手序号');
zlabel('速度/(m·s^{-1})');

v5_head = 2/max(max(v_mat))*v;
disp('问题5龙头速度: ');
disp(vpa(floor(v5_head*1e6)*1e-6,8));
