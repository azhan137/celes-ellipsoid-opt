%%
%for computing adjoint
%dFoM/dR_j = dFoM/dE_k*dE_k/db_i,n*db_i,n/dR_j
%adjointDerivative is 
%dFoM/dE_k*dE_k/db_i,n = (abs(E1)^2+abs(E2)^2)*conj(E_k)*N_i,n

function adjointDerivative = compute_adjoint_chromatic_intensity_derivative(simulation,E1,E2,points)

conjE1 = conj(E2).';
conjE2 = conj(E1).';
intensityE = abs(conjE1).^2 + abs(conjE2).^2;
%conjE = 2*real(conjE);
N = compute_bessel_value(simulation,points);

adjointDerivative = intensityE.*conjE1.*N;
adjointDerivative = squeeze(sum(squeeze(sum(adjointDerivative,1)),1));