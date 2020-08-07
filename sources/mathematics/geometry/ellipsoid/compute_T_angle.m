%%compute T-matrix of an ellipsoid
%ellip_params: 3x1 matrix: a,b,c along x,y,z
%ni: refractive index of surrounding medium
%ns: refractive index of particle
%lambda: incident wavelength
%Nphi, Ntheta: number of degrees for quadrature for phi, theta respectively
%lmax: maximum multipole order

%compute T,dT in lab frame

function [T,dT] = compute_T_angle(lmax,Ntheta,Nphi,ellip_params,ni,ns,lambda)
a = ellip_params(1,1);
b = ellip_params(1,2);
c = ellip_params(1,3);
phi_rot = ellip_params(1,4);

[Q,dQ] = compute_Q_angle(lmax,Ntheta,Nphi,a,b,c,phi_rot,ni,ns,lambda,3);
[rQ,drQ] = compute_Q_angle(lmax,Ntheta,Nphi,a,b,c,phi_rot,ni,ns,lambda,1);
Qinv = inv(Q);
T = rQ*Qinv;
dT = zeros([size(T),4]);
dT(:,:,1) = (drQ(:,:,1)-T*dQ(:,:,1))*Qinv;
dT(:,:,2) = (drQ(:,:,2)-T*dQ(:,:,2))*Qinv;
dT(:,:,3) = (drQ(:,:,3)-T*dQ(:,:,3))*Qinv;
dT(:,:,4) = (drQ(:,:,4)-T*dQ(:,:,4))*Qinv;