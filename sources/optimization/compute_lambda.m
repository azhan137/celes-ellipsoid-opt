%compute Lambda or adjointFields
%follows from
%L^T = adjointDerivative*inv(M)
%code calculates the following
%L = inv(M^T)*adjointDerivative^T

function adjointFields = compute_lambda(simulation,points)

adjointDerivative = gather(compute_adjoint_derivative(simulation,points).');
adjointDerivative = adjointDerivative(:);
masterMatrixTranspose = prepareM(simulation).';

% adjointFields = masterMatrixTranspose*adjointDerivative;
% % adjointFields = reshape(adjointFields,[simulation.numerics.nmax,simulation.input.particles.number]);
adjointFields = masterMatrixTranspose\adjointDerivative;
adjointFields = reshape(adjointFields,[simulation.numerics.nmax,simulation.input.particles.number]);