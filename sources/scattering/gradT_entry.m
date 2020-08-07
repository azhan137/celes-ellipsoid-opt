%======================================================================
%> @brief Evaluate the gradient of Mie coefficients with respect to radius
%> of the sphere from:
%> Grainger R., Lucas J., Thomas G., Ewen G., "Calculation of Mie
%> derivatives," Appl. Opt. 43, 5386-5393 (2004)
%>
%> @param       tau (int): SVWF polarization (1=TE, 2=TM)
%>
%> @param       l (int): SVWF degree (polar quantum number)
%>
%> @param       km (float): wavenumber in surrounding medium
%>
%> @param       kS (float): wavenumber inside particle
%>
%> @param       R (float): radius of sphere
%>
%> @retval      grad_Q(float): gradient of Mie coefficient

%> To do: internal field derivatives
%======================================================================
function grad_Q = gradT_entry(tau,l,kM,kS,R,varargin)

% The conventions are according to Bohren and Huffman's textbook.
% There is a twist: for tau=1 we have Q=-b and for tau=2 we have Q=-a in
% the case of scattered field Mie coefficients.
% This has to do with the definition of a and b in BH, see also
% Mishchenko: "Scattering, Absorption and Emission of Light by Small
% Paritcles", equations (5.42) and (5.43).

m = kS/kM;
x = kM*R;
mx = kS*R;

Amx = log_ricc_bessel_deriv(l,mx);
rhx = ricc_bessel(3,l,x);
rhx2 = ricc_bessel(3,l-1,x);

grad_a = 1i*kM*((Amx).^2*(1/m^2-1)+l*(l+1)*(1/mx^2-1/x^2))/((Amx/m+l/x)*rhx-rhx2)^2;
grad_b = 1i*kM*(1-m^2+l*(l+1)*((m/mx)^2-1/x^2))/((m*Amx+l/x)*rhx-rhx2)^2;

if(isempty(varargin))
    varargin={'scattered'};
end

switch varargin{1}
    case 'scattered'
        if tau==1
            grad_Q = grad_b;
        else
            grad_Q = grad_a;
        end
    otherwise
        error('only scattered implemented')
end
