function P = generate_dot(theta1,theta2,a)
% 通过板凳前后把手对应的角度求出板凳四个顶点的位置
    P = zeros(4,2);
    p1 = cal_xy(theta1,a);
    p2 = cal_xy(theta2,a);
    n1 = p1-p2;
    n0 = n1/norm(n1);
    p3 = p1+n0*0.275+n0*[0 1;-1 0]*0.15;
    p4 = p1+n0*0.275-n0*[0 1;-1 0]*0.15;
    p5 = p2-n0*0.275+n0*[0 1;-1 0]*0.15;
    p6 = p2-n0*0.275-n0*[0 1;-1 0]*0.15;
    P = [p3;p4;p5;p6];

end