function flag = judge_in(p,theta1,theta2,a)
% 判断一个点是否在某板凳的矩形区域内
    flag = 0;
    P = generate_dot(theta1,theta2,a);
    p1 = P(1,:);
    p2 = P(2,:);
    p3 = P(3,:);
    n0 = p-p1;
    n1 = p2-p1;
    n2 = p3-p1;
    if (dot(n0,n1)/dot(n1,n1)>0)&&(dot(n0,n1)/dot(n1,n1)<1)&&(dot(n0,n2)/dot(n2,n2)>0)&&(dot(n0,n2)/dot(n2,n2)<1)
        flag = 1;
    end

end