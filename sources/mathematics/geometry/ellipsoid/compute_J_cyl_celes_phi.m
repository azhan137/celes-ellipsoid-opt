%compute J matrices 11,12,21,22 using gaussian quadrature for computation
%   of full T matrix of an circular cylinder
%formulas are derived by Alan Zhan
%using piece-wise parameterization
%surface is parameterized as
%row(theta,phi) = h/2 1/cos(theta) for 0 to theta+
%row(theta,phi) = r/sin(theta) for theta+ to pi/2
%where theta+ is defined as arctan(2a/h)
%a is cylinder radius, h is cylinder height

%n_hat(theta,phi) = row_hat - tan(theta)theta_hat for 0 to theta+
%n_hat(theta,phi) = row_hat - cot(theta)theta_hat for theta+ to pi/2


%lmax: integer specifying maximum order
%Ntheta: number of points for Gaussian quadrature along polar
%Nphi: number of points for Gaussian quadrature along azimuthal
%a,h radius and height of cylinder
%ni,ns: refractive index of medium and particle respectively
%lambda: wavelength of light
%nu: internal (1) or external (3)
%

function [J11,J12,J21,J22,dJ11,dJ12,dJ21,dJ22] = compute_J_cyl_celes_phi(lmax,Ntheta,Nphi,a,h,n0,ns,lambda,nu)
%preallocate memory for J (nmax/2 x nmax/2)
nmax = jmult_max(1,lmax);
J11 = zeros(nmax/2);
J12 = zeros(nmax/2);
J21 = zeros(nmax/2);
J22 = zeros(nmax/2);
%preallocate for dJ (nmax/2 x nmax/2 x 2) 1:r, 2:h
dJ11 = zeros(nmax/2,nmax/2,2);
dJ12 = zeros(nmax/2,nmax/2,2);
dJ21 = zeros(nmax/2,nmax/2,2);
dJ22 = zeros(nmax/2,nmax/2,2);
%setup for gaussian quadrature in two dimensions
%first find angle where surface transitons
theta_edge = atan(2*a/h);
%then produce weights and integration points
[theta_circ,wtheta_circ] = generate_gauss_weights_abscissae(Ntheta,0,theta_edge);
[theta_body,wtheta_body] = generate_gauss_weights_abscissae(Ntheta,theta_edge,pi/2);
%combine the two
theta_i = [theta_body;theta_circ];
wtheta_i = [wtheta_body;wtheta_circ];
[phi_i,wphi_i] = generate_gauss_weights_abscissae(Nphi,0,pi);
[theta_i,wtheta_i] = generate_gauss_weights_abscissae(2*Ntheta,0,pi/2);
%generate angle and weight maps
[theta_map,phi_map] = meshgrid(theta_i,phi_i);
[wtheta_map,wphi_map] = meshgrid(wtheta_i,wphi_i);

weight_map = wtheta_map.*wphi_map;

%create indices for circle points and body points
[~,edge_index] = min(abs(theta_i-theta_edge));
body_index = 1:edge_index;
circ_index = (edge_index):2*Ntheta;

%k vector freespace and in scatterer
k = 2*pi*n0/lambda;
ks = 2*pi*ns/lambda;
%associated legendre polynomials
P_lm = legendre_normalized_angular(theta_map,lmax);
[Pi_lm,Tau_lm] = spherical_functions_angular(theta_map,lmax);
ct = cos(theta_map);
st = sin(theta_map);
tan_t = tan(theta_map(:,circ_index));
cot_t = cot(theta_map(:,body_index));

%parameterize the two radii:
r_circ = h/2./ct(:,circ_index);
r_body = a./st(:,body_index);
%combine into a single object
r = zeros(size(theta_map));
r(:,body_index) = r_body;
r(:,circ_index) = r_circ;
norm_factor_circ = 1;%./ct(:,circ_index);
norm_factor_body = 1;%./st(:,body_index);

for li = 1:lmax
    %precompute bessel functions for li,mi and derivative
    b_li = sph_bessel(nu,li,k*r);
    db_li = d1Z_Z_sph_bessel(nu,li,k*r);
    db2_li = (li+li.^2-(k*r).^2)./(k*r).*b_li;
    d1b_li = d1Z_sph_bessel(nu,li,k*r);
    for mi = -li:li
        %compute index ni
        ni = multi2single_index(1,1,li,mi,lmax);
        %get spherical harmonics (p,pi,tau) and derivatives
        p_limi = P_lm{li+1,abs(mi)+1};
        pi_limi = Pi_lm{li+1,abs(mi)+1};
        tau_limi = Tau_lm{li+1,abs(mi)+1};
        
        for lp = 1:lmax
            %precompute bessel functions for lp,mp, and derivative
            %for circular surface
            j_lp = sph_bessel(1,lp,ks*r);
            dj_lp = d1Z_Z_sph_bessel(1,lp,ks*r);
            dj2_lp = (lp+lp.^2-(ks*r).^2)./(ks*r).*j_lp;
            d1j_lp = d1Z_sph_bessel(1,lp,ks*r);
            lfactor = 1/sqrt(li*(li+li)*lp*(lp+1));
            for mp = -lp:lp
                %compute index np
                np = multi2single_index(1,1,lp,mp,lmax);
                %get spherical harmonics (p,pi,tau) and derivatives
                p_lpmp = P_lm{lp+1,abs(mp)+1};
                pi_lpmp = Pi_lm{lp+1,abs(mp)+1};
                tau_lpmp = Tau_lm{lp+1,abs(mp)+1};
              
                %setup selection rules for J11,J22
                selection_rules_1122 = selection_rules(li,mi,lp,mp,1);
                
                %setup selection rules for J12,J21
                selection_rules_1221 = selection_rules(li,mi,lp,mp,2);
                
                %phi exponential phase factor integrated into sinc
                phi_exp = exp(1i*(mp-mi)*phi_map);
                
                %compute J11,J22
                if selection_rules_1122 ~= 0
                    prefactor = selection_rules_1122*lfactor*phi_exp;
                    ang = mp*pi_lpmp.*tau_limi+mi*pi_limi.*tau_lpmp;
                    J11_r = -1i*weight_map.*prefactor.*r.^2.*st.*j_lp.*b_li.*ang;
                    J11(ni,np) = sum(sum(J11_r(:,circ_index)./norm_factor_circ))+sum(sum(J11_r(:,body_index)./norm_factor_body));
                    J22_r = -1i*prefactor.*weight_map.*st/k/ks.*dj_lp.*db_li.*ang;
                    J22_t = -1i*prefactor.*weight_map.*st/k/ks.*(lp*(lp+1)*mi*pi_limi.*p_lpmp.*j_lp.*db_li+li*(li+1)*mp.*pi_lpmp.*p_limi.*dj_lp.*b_li);
                    J22(ni,np) = sum(sum(J22_r(:,circ_index)./norm_factor_circ)+sum(J22_t(:,circ_index).*tan_t./norm_factor_circ))+sum(sum(J22_r(:,body_index)./norm_factor_body)+sum(J22_t(:,body_index).*-cot_t./norm_factor_body));            
                end
                %compute J12,J21
                if selection_rules_1221 ~= 0
                    prefactor = selection_rules_1221*lfactor*phi_exp;
                    J12_r = prefactor.*weight_map/k.*r.*st.*j_lp.*db_li.*(mi*mp*pi_limi.*pi_lpmp+tau_limi.*tau_lpmp);
                    J12_t = prefactor.*weight_map/k.*r.*st*li*(li+1).*j_lp.*b_li.*p_limi.*tau_lpmp;
                    J12(ni,np) = sum(sum(J12_r(:,circ_index)./norm_factor_circ)+sum(J12_t(:,circ_index).*tan_t./norm_factor_circ))+sum(sum(J12_r(:,body_index)./norm_factor_body)+sum(J12_t(:,body_index).*-cot_t./norm_factor_body));
                    J21_r = -prefactor.*weight_map/ks.*r.*st.*dj_lp.*b_li.*(mi*mp*pi_limi.*pi_lpmp+tau_limi.*tau_lpmp);
                    J21_t = -prefactor.*weight_map/ks.*r.*st*lp*(lp+1).*j_lp.*b_li.*p_lpmp.*tau_limi;
                    J21(ni,np) = sum(sum(J21_r(:,circ_index)./norm_factor_circ)+sum(J21_t(:,circ_index).*tan_t./norm_factor_circ))+sum(sum(J21_r(:,body_index)./norm_factor_body)+sum(J21_t(:,body_index).*-cot_t./norm_factor_body));
                end
            end
        end
    end
end
end
        

function selection_rules = selection_rules(l,m,lp,mp,diag_switch)
if diag_switch == 1
    selection_rules = (-1)^(m)*(1+(-1)^(mp-m))*(1+(-1)^(lp+l+1));
elseif diag_switch == 2
    selection_rules = (-1)^(m)*(1+(-1)^(mp-m))*(1+(-1)^(lp+l));
else
    disp('not supported, only valid inputs for diag_switch are 1,2');
end
end
                  
                