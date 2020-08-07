%======================================================================
%> @brief Riccati Bessel or Hankel function
%>
%> @param nu (int): selects between Riccati Bessel (nu=1) and Hankel (nu=3)
%> function
%> @param l (int): order of Riccati Bessel/Hankel function
%> @param Z (complex float array): argument 
%>
%> @retval rb (complex float array): Riccati Bessel or Hankel function
%======================================================================

function rb = ricc_bessel(nu,l,Z)
%function rb = ricc_bessel(nu,l,Z)
% riccati bessel function
% nu = 1: riccati bessel of the first kind (j) order l
% nu = 3: riccati hankel of the first kind (h1) order l
if nu == 1
    rb = Z*sph_bessel(nu,l,Z);
elseif nu == 3
    rb = -Z*sph_bessel(nu,l,Z);
else
    error('nu must be 1 or 3');
end
