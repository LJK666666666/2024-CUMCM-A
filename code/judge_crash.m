function flag = judge_crash(theta_list,a)
% 判断是否发生碰撞
    flag = 0;
    dot_mat = zeros(223,4,2);
    for i=1:223
        dot_mat(i,:,:) = generate_dot(theta_list(i),theta_list(i+1),a);
    end
    dot_list = [dot_mat(1,1,1),dot_mat(1,1,2);
                dot_mat(1,3,1),dot_mat(1,3,2);
                dot_mat(2,1,1),dot_mat(2,1,2)];
    for i=1:3
        for j=4:223
            dot = dot_list(i,:);
            
            if judge_in(dot,theta_list(j),theta_list(j+1),a)==1
                flag = 1;
            end
        end
    end

end