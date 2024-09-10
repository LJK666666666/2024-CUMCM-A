function n1_mat = cal_n1_all(theta_mat,flag_mat)
% 计算所有把手的切线的方向向量
    [x,y] = size(theta_mat);
    n1_mat = zeros(x,y,2);
    for i=1:x
        for j=1:y
            n1_mat(i,j,:) = cal_n1(theta_mat(i,j),flag_mat(i,j));
        end
    end

end