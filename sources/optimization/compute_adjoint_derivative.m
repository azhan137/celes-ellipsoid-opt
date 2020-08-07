%%
%for computing adjoint
%dFoM/dR_j = dFoM/dE_k*dE_k/db_i,n*db_i,n/dR_j
%adjointDerivative is dFoM/dE_k*dE_k/db_i,n = conj(E_k)*N_i,n

function adjointDerivative = compute_adjoint_derivative(simulation,E,points)

conjE = conj(E);
%conjE = 2*real(conjE);
N = compute_bessel_value(simulation,points);

adjointDerivative = squeeze(conjE.*N);
adjointDerivative = squeeze(sum(adjointDerivative,1));