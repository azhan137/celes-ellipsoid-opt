%compute full gradient (parallel)
%dI/dR_j = 2*real(L^T * d_j)
% d_j is RHS, L is adjoint coefficients

function adjointIntensityGradient = compute_parallel_adjoint_grad(simulation)
    %get particle information, number of particles, and indices
    numberParticles = simulation.input.particles.number;
    radiusArray = simulation.input.particles.radiusArrayIndex;
    radiusArray = radiusArray(:);
    %store required quantities
    Wb = obj.tables.coupledScatteringCoefficients;
    gradMie = obj.tables.gradMieCoefficients;
    initialField = obj.initialfieldCoefficients;
    adjointFields = (simulation.tables.adjointFields).';
    %allocate intensity gradient
    adjointIntensityGradient = zeros(numberParticles,1);

    parfor i = 1:numberParticles
        tempGradMie = zeros(size(tempGradMie));
        tempGradMie(i,:) = gradMie(radiusArray(i),:);
        RHS_i = tempGradMie.*Wb + initialField.*tempGradMie;
        adjointIntensityGradient(i) = sum(sum(adjointFields.*RHS_i));
    end
    
    adjointIntensityGradient = 2*real(adjointIntensityGradient);
end