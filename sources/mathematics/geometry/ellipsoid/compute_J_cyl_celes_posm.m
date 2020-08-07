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

function [J11,J12,J21,J22,dJ11,dJ12,dJ21,dJ22] = compute_J_cyl_celes_posm(lmax,Ntheta,a,h,ni,ns,lambda,nu)
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
%generate angle and weight maps
theta_map = theta_i;
weight_map = wtheta_i;

%create indices for circle points and body points
[~,edge_index] = min(abs(theta_map-theta_edge));
body_index = 1:edge_index;
circ_index = (edge_index+1):2*Ntheta;

%k vector freespace and in scatterer
k = 2*pi*ni/lambda;
ks = 2*pi*ns/lambda;
%associated legendre polynomials
P_lm = legendre_normalized_angular(theta_map,lmax);
[Pi_lm,Tau_lm] = spherical_functions_angular(theta_map,lmax);
ct = cos(theta_map);
st = sin(theta_map);
tan_t = tan(theta_map(circ_index));
cot_t = cot(theta_map(body_index));

%parameterize the two radii:
r_circ = h/2./ct(circ_index);
r_body = a./st(body_index);

%derivatives
dr_circ = 1/2./ct(circ_index);
dr_body = 1./st(body_index);
%combine into a single object
r = [r_body;r_circ];
norm_factor_circ = 1;%sqrt(1+tan_t.^2);
norm_factor_body = 1;%sqrt(1+cot_t.^2);

edgeda = 2*h/(h^2+4*a^2);
edgedh = -2*a/(h^2+4*a^2);

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
            lfactor = 1/sqrt(li*(li+1)*lp*(lp+1));
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
                
                %phi exponential phase factor integrated
                if mi ~= mp
                    phi_exp = -1i*(exp(1i*(mp-mi)*pi)-1)/(mp-mi);
                else
                    phi_exp = pi;
                end
                
                %compute J11,J22
                if selection_rules_1122 ~= 0
                    prefactor = selection_rules_1122*lfactor*phi_exp;
                    ang = mp*pi_lpmp.*tau_limi+mi*pi_limi.*tau_lpmp;
                    J11_r = -1i*weight_map.*prefactor.*r.^2.*st.*j_lp.*b_li.*ang;
                    dJ11da = 2*r_body.*j_lp(body_index).*b_li(body_index)+r_body.^2.*(ks*d1j_lp(body_index).*b_li(body_index)+k*d1b_li(body_index).*j_lp(body_index));
                    dJ11dh = 2*r_circ.*j_lp(circ_index).*b_li(circ_index)+r_circ.^2.*(ks*d1j_lp(circ_index).*b_li(circ_index)+k*d1b_li(circ_index).*j_lp(circ_index));
                    dJ11(ni,np,1) = sum(-1i*prefactor*weight_map(body_index).*st(body_index).*dJ11da.*dr_body.*ang(body_index));
                    dJ11(ni,np,2) = sum(-1i*prefactor*weight_map(circ_index).*st(circ_index).*dJ11dh.*dr_circ.*ang(circ_index));
                    J11(ni,np) = sum(J11_r(circ_index)./norm_factor_circ)+sum(J11_r(body_index)./norm_factor_body);
                    J22_r = -1i*prefactor.*weight_map.*st/k/ks.*dj_lp.*db_li.*ang;
                    J22db = lp*(lp+1)*mi*pi_limi.*p_lpmp;
                    J22dj = li*(li+1)*mp*pi_lpmp.*p_limi;
                    J22_t = -1i*prefactor.*weight_map.*st/k/ks.*(J22db.*j_lp.*db_li+J22dj.*dj_lp.*b_li);
                    dJ22edge = -1i/k/ks*st(edge_index)*(J22db(edge_index)*j_lp(edge_index)*db_li(edge_index)+J22dj(edge_index)*b_li(edge_index)*dj_lp(edge_index))*(st(edge_index)/ct(edge_index)+ct(edge_index)/st(edge_index));
                    dJ22da1 = -1i/k/ks*(ks*dj2_lp(body_index).*db_li(body_index)+k*db2_li(body_index).*dj_lp(body_index)).*dr_body.*st(body_index).*ang(body_index);
                    dJ22da2 = 1i/k/ks*cot_t.*st(body_index).*dr_body.*(J22db(body_index).*(ks*d1j_lp(body_index).*db_li(body_index)+k*j_lp(body_index).*db2_li(body_index))+J22dj(body_index).*(k*d1b_li(body_index).*dj_lp(body_index)+ks*b_li(body_index).*dj2_lp(body_index)));
                    dJ22dh1 = -1i/k/ks*(ks*dj2_lp(circ_index).*db_li(circ_index)+k*db2_li(circ_index).*dj_lp(circ_index)).*dr_circ.*st(circ_index).*ang(circ_index);
                    dJ22dh2 = -1i/k/ks*tan_t.*st(circ_index).*dr_circ.*(J22db(circ_index).*(ks*d1j_lp(circ_index).*db_li(circ_index)+k*j_lp(circ_index).*db2_li(circ_index))+J22dj(circ_index).*(k*d1b_li(circ_index).*dj_lp(circ_index)+ks*b_li(circ_index).*dj2_lp(circ_index)));
                    dJ22(ni,np,1) = sum(prefactor.*weight_map(body_index).*dJ22da1) + sum(prefactor.*weight_map(body_index).*dJ22da2) + prefactor*dJ22edge*edgeda;
                    dJ22(ni,np,2) = sum(prefactor.*weight_map(circ_index).*dJ22dh1) + sum(prefactor.*weight_map(circ_index).*dJ22dh2) + prefactor*dJ22edge*edgedh;
                    J22(ni,np) = sum(J22_r(circ_index)./norm_factor_circ)+sum(J22_r(body_index)./norm_factor_body)+sum(J22_t(circ_index).*tan_t./norm_factor_circ)+sum(J22_t(body_index).*-cot_t./norm_factor_body);            
                end
                %compute J12,J21
                if selection_rules_1221 ~= 0
                    prefactor = selection_rules_1221*lfactor*phi_exp;
                    ang = mi*mp.*pi_limi.*pi_lpmp+tau_limi.*tau_lpmp;
                    J12_r = prefactor.*weight_map/k.*r.*st.*j_lp.*db_li.*ang;
                    J12_t = prefactor.*weight_map/k.*r.*st*li*(li+1).*j_lp.*b_li.*p_limi.*tau_lpmp;
                    J12(ni,np) = sum(J12_r(circ_index)./norm_factor_circ)+sum(J12_r(body_index)./norm_factor_body)+sum(J12_t(circ_index).*tan_t./norm_factor_circ)+sum(J12_t(body_index).*-cot_t./norm_factor_body);
                    dJ12edge = li*(li+1)/k*r(edge_index)*st(edge_index)*j_lp(edge_index)*b_li(edge_index)*tau_lpmp(edge_index)*p_limi(edge_index)*(st(edge_index)/ct(edge_index)+ct(edge_index)/st(edge_index));
                    dJ12da1 = dr_body/k.*(j_lp(body_index).*db_li(body_index)+r_body.*(ks*d1j_lp(body_index).*db_li(body_index)+k*j_lp(body_index).*db2_li(body_index))).*st(body_index).*ang(body_index);
                    dJ12da2 = -li*(li+1)/k*dr_body.*(j_lp(body_index).*b_li(body_index)+r_body.*(ks*d1j_lp(body_index).*b_li(body_index)+k*j_lp(body_index).*d1b_li(body_index))).*cot_t.*st(body_index).*tau_lpmp(body_index).*p_limi(body_index);
                    dJ12dh1 = dr_circ/k.*(j_lp(circ_index).*db_li(circ_index)+r_circ.*(ks*d1j_lp(circ_index).*db_li(circ_index)+k*j_lp(circ_index).*db2_li(circ_index))).*st(circ_index).*ang(circ_index);
                    dJ12dh2 = li*(li+1)/k*dr_circ.*(j_lp(circ_index).*b_li(circ_index)+r_circ.*(ks*d1j_lp(circ_index).*b_li(circ_index)+k*j_lp(circ_index).*d1b_li(circ_index))).*tan_t.*st(circ_index).*tau_lpmp(circ_index).*p_limi(circ_index);
                    dJ12(ni,np,1) = sum(prefactor.*weight_map(body_index).*dJ12da1)+sum(prefactor.*weight_map(body_index).*dJ12da2)+prefactor*dJ12edge*edgeda;
                    dJ12(ni,np,2) = sum(prefactor.*weight_map(circ_index).*dJ12dh1)+sum(prefactor.*weight_map(circ_index).*dJ12dh2)+prefactor*dJ12edge*edgedh;   
                    J21_r = -prefactor.*weight_map/ks.*r.*st.*dj_lp.*b_li.*ang;
                    J21_t = -prefactor.*weight_map/ks.*r.*st*lp*(lp+1).*j_lp.*b_li.*p_lpmp.*tau_limi;
                    dJ21edge = -lp*(lp+1)/ks*r(edge_index)*st(edge_index)*j_lp(edge_index)*b_li(edge_index)*tau_limi(edge_index)*p_lpmp(edge_index)*(st(edge_index)/ct(edge_index)+ct(edge_index)/st(edge_index));
                    dJ21da1 = -dr_body/ks.*(b_li(body_index).*dj_lp(body_index)+r_body.*(k*d1b_li(body_index).*dj_lp(body_index)+ks*dj2_lp(body_index).*b_li(body_index))).*st(body_index).*ang(body_index);
                    dJ21da2 = lp*(lp+1)/ks*dr_body.*(j_lp(body_index).*b_li(body_index)+r_body.*(ks*d1j_lp(body_index).*b_li(body_index)+k*d1b_li(body_index).*j_lp(body_index))).*cot_t.*st(body_index).*tau_limi(body_index).*p_lpmp(body_index);
                    dJ21dh1 = -dr_circ/ks.*(b_li(circ_index).*dj_lp(circ_index)+r_circ.*(k*d1b_li(circ_index).*dj_lp(circ_index)+ks*dj2_lp(circ_index).*b_li(circ_index))).*st(circ_index).*ang(circ_index);
                    dJ21dh2 = -lp*(lp+1)/ks*dr_circ.*(j_lp(circ_index).*b_li(circ_index)+r_circ.*(ks*d1j_lp(circ_index).*b_li(circ_index)+k*d1b_li(circ_index).*j_lp(circ_index))).*tan_t.*st(circ_index).*tau_limi(circ_index).*p_lpmp(circ_index);
                    dJ21(ni,np,1) = sum(prefactor.*weight_map(body_index).*dJ21da1)+sum(prefactor.*weight_map(body_index).*dJ21da2)+prefactor*dJ21edge*edgeda;
                    dJ21(ni,np,2) = sum(prefactor.*weight_map(circ_index).*dJ21dh1)+sum(prefactor.*weight_map(circ_index).*dJ21dh2)+prefactor*dJ21edge*edgedh;
                    J21(ni,np) = sum(J21_r(circ_index)./norm_factor_circ)+sum(J21_r(body_index)./norm_factor_body)+sum(J21_t(circ_index).*tan_t./norm_factor_circ)+sum(J21_t(body_index).*-cot_t./norm_factor_body);
                    if ni == 10
                        abc = 2;
                    end
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
                    
                
                