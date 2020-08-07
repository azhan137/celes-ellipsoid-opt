function [T,dT] = compute_T_p(Q, rQ, dQ, drQ)

Qinv = inv(Q);

T_particle = rQ*Qinv;
dT_particle  = zeros([size(T_particle),2]);

dT_particle(:,:,1) = (drQ(:,:,1)-T_particle*dQ(:,:,1))*Qinv;
dT_particle(:,:,2) = (drQ(:,:,2)-T_particle*dQ(:,:,2))*Qinv;

T = T_particle;
dT = dT_particle;