function theta_n1 = cal_theta_n1_reverse(theta_n,d,a)
% 计算盘出螺线上下一个把手对应的角度
    theta_n1 = theta_n;
    lambda = pi/8;
    epsilon = 1e-10;
    while lambda>epsilon
        theta_n1 = theta_n1-lambda;
        if theta_n^2+theta_n1^2-2*theta_n*theta_n1*cos(theta_n1-theta_n)>d^2/a^2
            theta_n1 = theta_n1+lambda;
            lambda = lambda/2;
        end
    end            
end