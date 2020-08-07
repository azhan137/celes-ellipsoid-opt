radius = [200,200.1];
Ntheta = 30;
lmax = 4;
ni = 1;
ns = 2;
wavelength = 1000;

for i = 1:length(radius)
    [T,dT] = compute_T_cyl(lmax,Ntheta,[radius(i),460],ni,ns,wavelength);
    figure
    imagesc(abs(T))
    colorbar
    title(num2str(radius(i)))
    figure
    imagesc(abs(dT(:,:,1)))
    colorbar
    figure
    imagesc(abs(dT(:,:,2)))
    colorbar
    
end

lambda = wavelength;

[J1p,J2p,J3p,J4p] = compute_J_cyl_celes_phi(lmax,Ntheta,Ntheta,300,800,ni,ns,lambda,3);
[J11,J21,J31,J41] = compute_J_cyl_celes_posm(lmax,Ntheta,300,800,ni,ns,lambda,3);

figure
subplot(2,2,1)
imagesc(abs(J1p-J11))
colorbar
subplot(2,2,2)
imagesc(abs(J2p-J21))
colorbar
subplot(2,2,3)
imagesc(abs(J3p-J31))
colorbar
subplot(2,2,4)
imagesc(abs(J4p-J41))
colorbar

% [T_ellip, ~] = compute_T(lmax,Ntheta,Ntheta,[600,600,300],1,1,1000);

% [J1,J2,J3,J4] = compute_J_cyl_celes_posm(lmax,Ntheta,600,300,1,1,1000,3);
% [Ja,Jb,Jc,Jd] = compute_J_cyl_celes_unnorm(lmax,Ntheta,600,300,1,1,1000,3);
% [J1p,J2p,J3p,J4p] = compute_J_cyl_celes_phi(lmax,Ntheta,Nphi,600,300,1,1,1000,3);
% [J1e,J2e,J3e,J4e] = compute_J_ellip_celes_posm(lmax,Ntheta,Ntheta,300,300,600,1,1,1000,3);
% [rJ1e,rJ2e,rJ3e,rJ4e] = compute_J_ellip_celes_posm(lmax,Ntheta,Ntheta,300,300,600,1,1,1000,1);
% [Qe,~] = compute_Q(lmax,Ntheta,Ntheta,300,300,600,1,1,1000,3);
% [rQe,~] = compute_Q(lmax,Ntheta,Ntheta,300,300,600,1,1,1000,1);