function P = generate_dot_(p1,p2)
% 通过板凳前后把手的位置求出板凳四个顶点的位置

    P = zeros(4,2);
    p1 = reshape(p1,1,2,1);
    p2 = reshape(p2,1,2,1);
    n1 = p1-p2;
    n0 = n1/norm(n1);
    p3 = p1+n0*0.275+n0*[0 1;-1 0]*0.15;
    p4 = p1+n0*0.275-n0*[0 1;-1 0]*0.15;
    p5 = p2-n0*0.275+n0*[0 1;-1 0]*0.15;
    p6 = p2-n0*0.275-n0*[0 1;-1 0]*0.15;
    P = [p3;p4;p5;p6];

end