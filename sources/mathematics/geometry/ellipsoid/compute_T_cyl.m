function [T,dT] = compute_T_cyl(lmax,Ntheta,cyl_params,ni,ns,lambda)
a = cyl_params(1,1);
h = cyl_params(1,2);

[Q,dQ] = compute_Q_cyl(lmax,Ntheta,a,h,ni,ns,lambda,3);
[rQ,drQ] = compute_Q_cyl(lmax,Ntheta,a,h,ni,ns,lambda,1);

Qinv = inv(Q);

T_particle = rQ*Qinv;
dT_particle  = zeros([size(T_particle),2]);

dT_particle(:,:,1) = (drQ(:,:,1)-T_particle*dQ(:,:,1))*Qinv;
dT_particle(:,:,2) = (drQ(:,:,2)-T_particle*dQ(:,:,2))*Qinv;

T = T_particle;
dT = dT_particle;