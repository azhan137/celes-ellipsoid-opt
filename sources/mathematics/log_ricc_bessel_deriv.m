%======================================================================
%> @brief Logarithmic derivative of Riccati Bessel function of the first
%> kind
%> reference: Grainger R., Lucas J., Thomas G., Ewen G., "Calculation of Mie
%> derivatives," Appl. Opt. 43, 5386-5393 (2004)
%> @param l (int): order of Riccati Bessel function
%> @param Z (complex float array): argument 
%>
%> @retval An (complex float array): Logarithmic derivative of Riccati Bessel
%======================================================================

function An = log_ricc_bessel_deriv(l,Z)
%function An = log_ricc_bessel_deriv(l,Z)
% calculates logarithmic derivative used in derivatives of mie coefficients
% An = ricc_bessel(l-1,Z)-l/Z*ricc_bessel(l,z)/(ricc_bessel(l,z))
An = (ricc_bessel(1,l-1,Z)-l/Z*ricc_bessel(1,l,Z))/ricc_bessel(1,l,Z);
