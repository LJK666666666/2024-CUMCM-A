function n1 = cal_n1(theta,flag)
% 计算把手的切线的方向向量
    if flag == 1 | flag == 4
        n1 = [cos(theta)-theta.*sin(theta),sin(theta)+theta.*cos(theta)];
    else
        n1 = [sin(theta),-cos(theta)];
    end

end