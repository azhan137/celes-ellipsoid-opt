%compute full gradient
%dI/dR_j = 2*real(L^T * d_j)
%where d_j is the rhs of the equation M*db_i,n/dR_j = rhs

function adjointIntensityGradient = compute_adjoint_chromatic_grad_intensity(simulation)

numberParticles = simulation.input.particles.number;
radiusArray = simulation.input.particles.radiusArrayIndex;
radiusArray = radiusArray(:);

adjointFields = simulation.tables.adjointFields;
adjointFields = adjointFields.';

adjointIntensityGradient = zeros(numberParticles,1);

for i = 1:numberParticles
    rhs_i = gather(simulation.rightHandSideOpt(i,radiusArray(i)));
    adjointIntensityGradient(i) = sum(sum(adjointFields.*rhs_i));
end

adjointIntensityGradient = 2*real(adjointIntensityGradient);
    
    
