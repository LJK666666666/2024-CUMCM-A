function theta_n_list = cal_theta_all(theta_1,a)
% 计算所有把手对应的角度
    theta_n_list = [theta_1];
    
    for i=1:223
        theta_n1 = cal_theta_n1(theta_n_list(end),cal_d(i),a);
        theta_n_list = [theta_n_list;theta_n1];
    end
end