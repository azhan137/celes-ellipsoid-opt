r_0 = 200;
h_0 = 460;
lambda = 1000;
ns = 2;
ni = 1;

step_size = 0.01;
Ntheta = 100;
lmax = 4;

[J11, J12, J21, J22, dJ11, dJ12, dJ21, dJ22] = compute_J_cyl_celes_posm(lmax,Ntheta,200,460,ni,ns,lambda,3);
[J11j, J12j, J21j, J22j, dJ11j, dJ12j, dJ21j, dJ22j] = compute_J_cyl_celes_posm(lmax,Ntheta,200,460,ni,ns,lambda,1);

dJ11a = dJ11(:,:,1);
dJ11h = dJ11(:,:,2);
dJ12a = dJ12(:,:,1);
dJ12h = dJ12(:,:,2);
dJ21a = dJ21(:,:,1);
dJ21h = dJ21(:,:,2);
dJ22a = dJ22(:,:,1);
dJ22h = dJ22(:,:,2);

[T_0,dT_0] = compute_T_cyl(lmax,Ntheta,[r_0,h_0],ni,ns,lambda);
[T_r,dT_r] = compute_T_cyl(lmax,Ntheta,[r_0+step_size,h_0],ni,ns,lambda);
[T_h,dT_h] = compute_T_cyl(lmax,Ntheta,[r_0,h_0+step_size],ni,ns,lambda);

numerical_h = (T_h-T_0)/step_size;
numerical_r = (T_r-T_0)/step_size;

figure
subplot(3,1,1)
imagesc(abs(numerical_r))
colorbar
axis square
subplot(3,1,2)
imagesc(abs(dT_0(:,:,1)))
colorbar
axis square
subplot(3,1,3)
imagesc(abs(dT_0(:,:,1)-numerical_r))
colorbar
axis square

figure
subplot(3,1,1)
imagesc(abs(numerical_h))
colorbar
axis square
subplot(3,1,2)
imagesc(abs(dT_0(:,:,2)))
colorbar
axis square
subplot(3,1,3)
imagesc(abs(dT_0(:,:,2)-numerical_h))
colorbar
axis square

r_i = r_0;
T_u = T_0;
% for i = 1:50
%     [T_i,dT_i] = compute_T_cyl(lmax,Ntheta,[r_i,h_0],ni,ns,lambda);
%     r_i = r_i + step_size;
%     T_u = T_u + dT_i(:,:,1)*step_size;
% end
% 
% figure
% imagesc(abs(T_u-T_i))
% colorbar