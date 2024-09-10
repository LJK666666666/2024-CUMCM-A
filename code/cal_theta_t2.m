function theta_t2 = cal_theta_t2(theta_t1,v,t,a,direction)
% 计算把手经过时间t后对应的角度
    
    lambda = theta_t1/4;
    epsilon = 1e-12;
    % epsilon = 1e-10;
    theta_t2 = theta_t1;
    if direction==1
        obj = cal_integral(theta_t1)-v.*t./a.*2;
        while lambda>epsilon
            theta_t2 = theta_t2-lambda;
            if cal_integral(theta_t2)<obj
                theta_t2 = theta_t2+lambda;
                lambda = lambda/2;
            end
        end
    elseif direction==-1
        obj = cal_integral(theta_t1)+v.*t./a.*2;
        while lambda>epsilon
            theta_t2 = theta_t2+lambda;
            if cal_integral(theta_t2)>obj
                theta_t2 = theta_t2-lambda;
                lambda = lambda/2;
            end
        end
    end
end
