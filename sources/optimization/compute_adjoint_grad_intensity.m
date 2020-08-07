%compute full gradient
%dI/dR_j = 2*real(L^T * d_j)
%where d_j is the rhs of the equation M*db_i,n/dR_j = rhs

%%varargin: 
%%empty for spheres
%%for ellipsoids (1,2,3,4) -> (a,b,c,phi) gradient
%%for cylinders (1,2) -> (r,h) gradient


function adjointIntensityGradient = compute_adjoint_grad_intensity(simulation,varargin)

numberParticles = simulation.input.particles.number;

adjointFields = simulation.tables.adjointFields;
adjointFields = adjointFields.';

switch simulation.input.particles.type
    case 'sphere'
        radiusArray = simulation.input.particles.singleParticleArrayIndex;
        radiusArray = radiusArray(:);
        
        adjointIntensityGradient = zeros(numberParticles,1);
        
        for i = 1:numberParticles
            rhs_i = gather(simulation.rightHandSideOpt(i,radiusArray(i)));
            adjointIntensityGradient(i) = sum(sum(adjointFields.*rhs_i));
        end
    case 'ellipsoid'
        particleArray = simulation.input.particles.singleParticleArrayIndex;
        particleArray = particleArray(:);
        
        adjointIntensityGradient = zeros(numberParticles,1);
        ellipsoidFlag = varargin{1};
        
        for i = 1:numberParticles
            rhs_i = gather(simulation.rightHandSideOpt(i,particleArray(i),ellipsoidFlag));
            adjointIntensityGradient(i) = sum(sum(adjointFields.*rhs_i));
        end
    case 'cylinder'
        particleArray = simulation.input.particles.singleParticleArrayIndex;
        particleArray = particleArray(:);
        
        adjointIntensityGradient = zeros(numberParticles,1);
        cylinderFlag = varargin{1};
        
        for i = 1:numberParticles
            rhs_i = gather(simulation.rightHandSideOpt(i,particleArray(i),cylinderFlag));
            adjointIntensityGradient(i) = sum(sum(adjointFields.*rhs_i));
        end
end
        

adjointIntensityGradient = 2*real(adjointIntensityGradient);
    
    
