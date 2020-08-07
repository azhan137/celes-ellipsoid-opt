%%test agreement between axial_rotation(T) and compute_T_angle where the
%%following:
%axial_rotation(T): compute T in particle frame, rotate to lab frame
%compute_T_angle: compute T in lab frame

lmax = 4;
Ntheta = 20;
Nphi = 20;
ellip_params = [1000,600,1000];
perp_ellip_params = [600,1000,1000];
ni = 1;
ns = 1.52;
lambda = 1550;
alpha = pi/2;
i = 0;
% 
% for alpha = 0:0.1:pi
% i = i+1;
%compute particle frame T and dT
[T_p,dT_p] = compute_T(lmax,Ntheta,Nphi,ellip_params,ni,ns,lambda);
%rotate from particle frame into lab frame
[T_rl,dT_rl] = axial_rotation(lmax,T_p,dT_p,-alpha);
%[T_pn,dT_pn] = axial_rotation(lmax,T_rl,dT_rl(:,:,1:3),-alpha);

%compute lab frame T and dT directly using different parameterization
[T_l,dT_l] = compute_T_angle(lmax,Ntheta,Nphi,[ellip_params,alpha],ni,ns,lambda);

%compute lab frame T and dT directly
[T_c,dT_c] = compute_T(lmax,Ntheta,Nphi,perp_ellip_params,ni,ns,lambda);
%get angular rotation matrix
[T_cf,dT_cf] = axial_rotation(lmax,T_c,dT_c,0);

figure(1)
imagesc(real(dT_rl(:,:,4)))
colorbar()
figure(2)
imagesc(real(dT_cf(:,:,4)))
colorbar()
% figure
% imagesc(abs(T_rl-T_l))
% colorbar
% title(num2str(alpha))
% angle_i = angle(T_rl(10,12)-T_l(10,12));
% angleList(i) = angle_i;
% max_i = abs(T_rl(10,12)-T_l(10,12));
% maxList(i) = max_i;
% real_i = real(T_rl(10,12)-T_l(10,12));
% realList(i) = real_i;
% imag_i = imag(T_rl(10,12)-T_l(10,12));
% imagList(i) = imag_i;
% end
% 
% figure
% plot(angleList)
% figure
% plot(maxList)