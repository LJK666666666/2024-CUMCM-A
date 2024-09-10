function p = cal_xy(theta_n_list,a)
% 通过螺线上把手对应的角度求出其位置
    x_list = a.*theta_n_list.*cos(theta_n_list);
    y_list = a.*theta_n_list.*sin(theta_n_list);
    p = [x_list,y_list];
end