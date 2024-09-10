function v2 = cal_v2(v1,n0,n1,n2)
% 通过上一个把手的速度推出下一个把手的速度

    n0 = reshape(n0,2,1,1);
    n1 = reshape(n1,2,1,1);
    n2 = reshape(n2,2,1,1);

    cos_yamma1 = abs(dot(n0,n1)/norm(n0)/norm(n1));
    cos_yamma2 = abs(dot(n0,n2)/norm(n0)/norm(n2));

    v2 = v1*cos_yamma1/cos_yamma2;

end