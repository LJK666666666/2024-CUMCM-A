function v_mat = cal_v_all(theta_mat,v)
% 计算所有把手所有时刻的速度

    num_t = size(theta_mat,2);

    tl = zeros(224,num_t,2);
    tl(:,:,1) = cos(theta_mat)-theta_mat.*sin(theta_mat);
    tl(:,:,2) = sin(theta_mat)+theta_mat.*cos(theta_mat);
    
    cl = zeros(223,num_t,2);
    cl(:,:,1) = theta_mat(1:223,:).*cos(theta_mat(1:223,:))-theta_mat(2:224,:).*cos(theta_mat(2:224,:));
    cl(:,:,2) = theta_mat(1:223,:).*sin(theta_mat(1:223,:))-theta_mat(2:224,:).*sin(theta_mat(2:224,:));
    
    v_mat = zeros(224,num_t);
    v_mat(1,:) = v;
    for j=1:num_t
        for i=2:224
            n1 = [cl(i-1,j,1),cl(i-1,j,2)];
            n2 = [tl(i-1,j,1),tl(i-1,j,2)];
            n3 = [tl(i,j,1),tl(i,j,2)];
            cos_yamma1 = abs(dot(n1,n2)/norm(n1)/norm(n2));
            cos_yamma2 = abs(dot(n1,n3)/norm(n1)/norm(n3));
            v_mat(i,j) = v_mat(i-1,j)*cos_yamma1/cos_yamma2;
        end
    end
end