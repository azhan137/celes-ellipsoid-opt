function grad_E = compute_electric_field_gradient(simulation,points)

Rlast = bsxfun(@plus,points,-simulation.input.particles.positionArray(simulation.input.particles.number,:));
rlast = sqrt(Rlast(:,1).^2+Rlast(:,2).^2+Rlast(:,3).^2);
rmax = max(rlast);  % estimate for maximal distance
resol = 1; % resolution of lookup
ri = 0:resol:(rmax+resol);
hi = cell(simulation.numerics.lmax,1);
dhi = cell(simulation.numerics.lmax,1);
for l=1:simulation.numerics.lmax
    hi{l} = simulation.numerics.deviceArray(single(sph_bessel(3,l,simulation.input.k_medium * ri)));
    dhi{l} = simulation.numerics.deviceArray(single(dx_xz(3,l,simulation.input.k_medium * ri)));
end
ri = simulation.numerics.deviceArray(single(ri));

grad_E = simulation.numerics.deviceArray(zeros([length(points),simulation.input.particles.number],'single'));

for grad_index = 1:simulation.input.particles.number
    for particle_index = 1:simulation.input.particles.number
        kM = simulation.input.k_medium;
        
        % relative positions
        R = simulation.numerics.deviceArray(single(bsxfun(@plus,points,-simulation.input.particles.positionArray(particle_index,:))));
        
        r = sqrt(R(:,1).^2+R(:,2).^2+R(:,3).^2);
        theta = acos( min(max(R(:,3)./r,-1),1));
        phi = atan2(R(:,2),R(:,1));
        
        e_r = R./[r,r,r];
        e_theta = [cos(theta).*cos(phi),cos(theta).*sin(phi),-sin(theta)];
        e_phi = [-sin(phi),cos(phi),r-r];
        
        kr = kM*r;
        % check if bessel lookup is sufficient
        if max(r)>ri(end)
            rnew = (ri(end)+resol):resol:((max(r)+resol)*1.4);
            for l=1:simulation.numerics.lmax
                hi{l} = [hi{l},simulation.numerics.deviceArray(single(sph_bessel(3,l,simulation.input.k_medium * gather(rnew))))];
                dhi{l} = [dhi{l},simulation.numerics.deviceArray(single(dx_xz(3,l,simulation.input.k_medium * gather(rnew))))];
            end
            ri = simulation.numerics.deviceArray(single([ri,rnew]));
        end
        
        % spherical functions
        [p_all] = legendre_normalized_angular(theta,simulation.numerics.lmax);
        [pi_all,tau_all] = spherical_functions_angular(theta,simulation.numerics.lmax);
        
        for l=1:simulation.numerics.lmax
            z=interp1(ri,hi{l},r,'linear');
            dxxz = interp1(ri,dhi{l},r,'linear');
            
            for m=-l:l
                P_lm = p_all{l+1,abs(m)+1};
                pi_lm = pi_all{l+1,abs(m)+1};
                tau_lm = tau_all{l+1,abs(m)+1};
                eimphi=exp(1i*m*phi);
                for tau=1:2
                    n = multi2single_index(1,tau,l,m,simulation.numerics.lmax);
                    % SVWFs
                    if tau==2  %select M
                        fac1 = 1/sqrt(2*l*(l+1)) * z .* eimphi;
                        fac2 = bsxfun(@times,1i*m*pi_lm,e_theta) - bsxfun(@times,tau_lm,e_phi);
                        N =  bsxfun(@times,fac1,fac2);
                    else %select N
                        fac1 = 1/sqrt(2*l*(l+1)) * eimphi;
                        term1 = bsxfun(@times,l*(l+1)*z./kr.*P_lm,e_r);
                        term22 = bsxfun(@times,tau_lm,e_theta);
                        term23 = bsxfun(@times,1i*m*pi_lm,e_phi);
                        term2 = bsxfun(@times,dxxz./kr,term22+term23);
                        N = bsxfun(@times,fac1,term1+term2);
                    end
                    grad_E(:,grad_index) = grad_E(:,grad_index) + squeeze(simulation.tables.gradScatteredFieldCoefficients(particle_index,n,grad_index)*N)';
                end
            end
        end
    end
end
grad_E = squeeze(grad_E);
        
        