clear all;
close all;

addpath(genpath('.'))

%%optimize a small array (10x10) of ellipsoids, changing their geometric 
%properties (axes a, b, c). Angle is not changed here, as it does not
%affect the performance in a good way currently. There is probably some
%problem with the way the gradient behaves when angular inputs are enabled.


%%%forward simulation only

simulation = celes_simulation;
particles = celes_particles2;
initialField = celes_initialField;
input = celes_input;
numerics = celes_numerics_v2;
solver = celes_solver;
inverseSolver = celes_solver;
preconditioner = celes_preconditioner_v2;
inversePreconditioner = celes_preconditioner_v2;

lmax = 3;
cuda_compile(lmax);
cuda_compile_T(lmax);

%particle properties
parameters = zeros(100,5);
a = ones(100,1)*125;
b = ones(100,1)*125;
c = ones(100,1)*400;
alpha = zeros(100,1);
%prevent overlap

max_rad = 180;
min_rad = 25;

%particle positions on a square grid
xpos = linspace(-10000,10000,10);
ypos = xpos';
[xx, yy] = meshgrid(xpos,ypos);
zz = zeros(length(a),1);
positions = zeros(length(a),3);
positions(:,1) = xx(:);
positions(:,2) = yy(:);
positions(:,3) = zz(:);
refractiveIndex = 2;

%total number of gradient descents to run
num_iterations = 50;
%max gradient step size
max_step_size = 100;
%optimize E-field intensity here
points = [0, 0, 20000];

parameters(:,1) = a;
parameters(:,2) = b;
parameters(:,3) = c;
parameters(:,4) = alpha;
parameters(:,5) = refractiveIndex*ones(size(a));

particles.parameterArray = parameters;
particles.positionArray = positions;
particles.type = 'ellipsoid';

%numeric properties
numerics.lmax = lmax;
numerics.particleDistanceResolution = 1;
numerics.gpuFlag = true;
numerics.polarAnglesArray = 0:pi/1e3:pi;
numerics.azimuthalAnglesArray = 0:pi/1e3:2*pi;

%solver properties
solver.type = 'BiCGStab';
solver.tolerance = 1e-3;
solver.maxIter = 1000;
solver.restart = 1000;
inverseSolver = solver;

%preconditioner properties
preconditioner.type = 'blockdiagonal';
inversePreconditioner.type = 'blockdiagonal';
numerics.partitionEdgeSizes = [5000,5000,6200];


%inptu into solver
solver.preconditioner = preconditioner;
inverseSolver.preconditioner = inversePreconditioner;

%initialFields
initialField.polarAngle = 0;
initialField.azimuthalAngle = 0;
initialField.polarization = 'TE';
initialField.beamWidth = 0;
initialField.focalPoint = [0,0,0];

%input properties minus wavelength
input.mediumRefractiveIndex = 1;

%put into simulation object;
input.initialField = initialField;
input.particles = particles;
simulation.input = input;
simulation.tables.pmax = simulation.input.particles.number;

numerics.solver = solver;
numerics.inverseSolver = inverseSolver;
simulation.numerics = numerics;
simulation.tables = celes_tables;
simulation.output = celes_output;
simulation.tables.nmax = simulation.numerics.nmax;
simulation.tables.pmax = simulation.input.particles.number;

simulation.input.wavelength = 650;

simulation = simulation.computeInitialFieldPower;
simulation = simulation.computeTranslationTable;
simulation.input.particles = simulation.input.particles.compute_maximal_particle_distance;

if strcmp(simulation.numerics.solver.preconditioner.type,'blockdiagonal')
    fprintf(1,'make particle partition ...');
    partitioning = make_particle_partion(simulation.input.particles.positionArray,simulation.numerics.partitionEdgeSizes);
    simulation.numerics.partitioning = partitioning;
    simulation.numerics.solver.preconditioner.partitioning = partitioning;
    simulation.numerics.inverseSolver.preconditioner.partitioning = partitioning;
    fprintf(1,' done\n');
    simulation = simulation.numerics.prepareW(simulation);
    simulation.numerics.solver.preconditioner.partitioningIdcs = simulation.numerics.partitioningIdcs;
    simulation.numerics.inverseSolver.preconditioner.partitioningIdcs = simulation.numerics.partitioningIdcs;
end

nest_params = parameters(:,1:2);

for i = 1:num_iterations
    
    simulation = simulation.computeParallelMieCoefficients;
    %Compute scattered field coefficients b_i,n
    simulation = simulation.numerics.solver.preconditioner.prepareM(simulation);

    simulation = simulation.computeInitialFieldCoefficients;
    simulation = simulation.computeScatteredFieldCoefficients();
    
    E1 = compute_scattered_field_opt(simulation, points);
    fom(i) = sum(abs(E1(:)).^2);
    
    fom(i)
    
    simulation = simulation.computeCoupledScatteringCoefficients;
    
    simulation = simulation.numerics.inverseSolver.preconditioner.prepareMt(simulation);
    simulation = simulation.computeAdjointCoefficients2(E1,points);
    
    %clear factored matrices
    simulation.numerics.inverseSolver.preconditioner.factorizedMasterMatrices = [];
    
    %compute gradient
    grad_a = compute_adjoint_grad_intensity(simulation,1);
    grad_b = compute_adjoint_grad_intensity(simulation,2);
    grad_c = compute_adjoint_grad_intensity(simulation,3);
    grad_na = grad_a/sqrt(sum(abs(grad_a(:).^2)));
    grad_nb = grad_b/sqrt(sum(abs(grad_b(:).^2)));
    grad_nc = grad_c/sqrt(sum(abs(grad_c(:).^2)));

    %normalize gradient
    full_grad = max_step_size*[grad_na, grad_nb, grad_nc];
    
    old_params = simulation.input.particles.parameterArray(:,1:3);
    new_params = old_params+full_grad;
    
    new_params(new_params > max_rad) = max_rad - 5;
    new_params(new_params < min_rad) = min_rad + 5;

    if mod(i,5) == 0
        save(filename, 'old_params', 'fom', 'points', 'positions', 'refractiveIndex');
    end
    
    simulation.input.particles.parameterArray(:,1:3) = new_params;    
end


[x,z] = meshgrid(-20000:20:20000,-1000:20:6000); y=x-x+000;
output.fieldPoints = [x(:),y(:),z(:)];
fieldPoints = [x(:),y(:),z(:)];
output.fieldPointsArrayDims = size(x);
simulation.output = output;
simulation=simulation.evaluateFields;
figure
plot_field(gca,simulation,'abs E','Total field',simulation.input.particles.parameterArray(:,1))
colormap(jet)
caxis([0 10])
colorbar
output.fieldPointsArrayDims = size(x);

points = [x(:), y(:), z(:)];
E = compute_scattered_field_opt(simulation, points)+compute_initial_field(simulation);
Ex = abs(E(:,1));
Ey = abs(E(:,2));

Ex = gather(reshape(Ex,121,121));
Ey = gather(reshape(Ey,121,121));

figure
imagesc(Ex)
colorbar
figure
imagesc(Ey)
colorbar

