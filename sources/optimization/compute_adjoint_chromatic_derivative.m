%%
%for computing adjoint
%dFoM/dR_j = dFoM/dE_k*dE_k/db_i,n*db_i,n/dR_j
%adjointDerivative is dFoM/dE_k*dE_k/db_i,n = conj(E_k)*N_i,n

function adjointDerivative = compute_adjoint_chromatic_derivative(simulation,E1,E2,points)

conjE1 = conj(E2).';
conjE2 = conj(E1).';
conjE = conjE1 + conjE2;
%conjE = 2*real(conjE);
N = compute_bessel_value(simulation,points);

adjointDerivative = conjE.*N;
adjointDerivative = squeeze(sum(squeeze(sum(adjointDerivative,1)),1));