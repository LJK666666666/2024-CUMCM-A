function integral = cal_integral(theta)
% 计算积分值

    integral = theta.*sqrt(theta.^2+1)+log(theta+sqrt(theta.^2+1));

end