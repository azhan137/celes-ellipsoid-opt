%%compute T-matrix of an ellipsoid
%ellip_params: 3x1 matrix: a,b,c along x,y,z
%ni: refractive index of surrounding medium
%ns: refractive index of particle
%lambda: incident wavelength
%Nphi, Ntheta: number of degrees for quadrature for phi, theta respectively
%lmax: maximum multipole order

%see pc waterman paper on extended boundary condition method

function [T,dT] = compute_T_q(lmax,Ntheta,Nphi,ellip_params,ni,ns,lambda)
a = ellip_params(1,1);
b = ellip_params(1,2);
c = ellip_params(1,3);
alpha = ellip_params(1,4);

[Q_p,dQ_p] = compute_Q(lmax,Ntheta,Nphi,a,b,c,ni,ns,lambda,3);
[rQ_p,drQ_p] = compute_Q(lmax,Ntheta,Nphi,a,b,c,ni,ns,lambda,1);

Q_l = zeros(size(Q_p));
dQ_l = zeros([size(Q_p),4]);
rQ_l = Q_l;
drQ_l = dQ_l;

for ti = 1:2
    for tp = 1:2
        for li = 1:lmax
            for lp = 1:lmax
                for mi = -li:li
                    ni = multi2single_index(1,ti,li,mi,lmax);
                    for mp = -lp:lp
                        np = multi2single_index(1,tp,lp,mp,lmax);
                        Q_l(ni,np) = exp(-1i*alpha*(mi-mp))*Q_p(ni,np);
                        dQ_l(ni,np,1:3) = exp(-1i*alpha*(mi-mp))*dQ_p(ni,np,:);
                        dQ_l(ni,np,4) = -1i*(mi-mp)*exp(-1i*alpha*(mi-mp))*Q_l(ni,np);
                        rQ_l(ni,np) = exp(-1i*alpha*(mi-mp))*rQ_p(ni,np);
                        drQ_l(ni,np,1:3) = exp(-1i*alpha*(mi-mp))*drQ_p(ni,np,:);
                        drQ_l(ni,np,4) = -1i*(mi-mp)*exp(-1i*alpha*(mi-mp))*rQ_l(ni,np);
                    end
                end
            end
        end
    end
end

dT = zeros(size(dQ_l));

Qinv_l = inv(Q_l);

T = rQ_l*Qinv_l;
dT(:,:,1) = (drQ_l(:,:,1)-T*dQ_l(:,:,1))*Qinv_l;
dT(:,:,2) = (drQ_l(:,:,2)-T*dQ_l(:,:,2))*Qinv_l;
dT(:,:,3) = (drQ_l(:,:,3)-T*dQ_l(:,:,3))*Qinv_l;
dT(:,:,4) = (drQ_l(:,:,4)-T*dQ_l(:,:,4))*Qinv_l;