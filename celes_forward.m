addpath(genpath('.'))

%%%forward simulation only

%% simulates 2x2 grid of spheroids and plots the electric fields a distance 
% away from them. Simple simulation for setting up a configuration in CELES

simulation = celes_simulation;
particles = celes_particles2;
initialField = celes_initialField;
input = celes_input;
numerics = celes_numerics_v2;
solver = celes_solver;
inverseSovler = celes_solver;
output = celes_output;
preconditioner = celes_preconditioner_v2;


lmax = 4;
cuda_compile(lmax);

%particle properties

parameters = zeros(4,4);
%%ellipsoid geometric parameters
a = ones(4,1)*150;
b = a;
c = a;
alpha = zeros(size(a));

xpos = linspace(-300,300,2);
ypos = xpos';
[xx, yy] = meshgrid(xpos,ypos);
zz = zeros(length(a),1);
positions = zeros(length(a),3);
positions(:,1) = xx(:);
positions(:,2) = yy(:);
positions(:,3) = zz(:);
refractiveIndex = 1.52;

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

%preconditioner properties
preconditioner.type = 'blockdiagonal';
numerics.partitionEdgeSizes = [15000,15000,6200];

%inptu into solver
solver.preconditioner = preconditioner;

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
simulation.numerics = numerics;
simulation.tables = celes_tables;
simulation.output = output;
simulation.tables.nmax = simulation.numerics.nmax;
simulation.tables.pmax = simulation.input.particles.number;

simulation.input.wavelength = 650;

simulation = simulation.computeInitialFieldPower;
simulation = simulation.computeTranslationTable;
simulation.input.particles = simulation.input.particles.compute_maximal_particle_distance;

%generate blockdiagonal preconditioner
if strcmp(simulation.numerics.solver.preconditioner.type,'blockdiagonal')
    fprintf(1,'make particle partition ...');
    partitioning = make_particle_partion(simulation.input.particles.positionArray,simulation.numerics.partitionEdgeSizes);
    simulation.numerics.partitioning = partitioning;
    simulation.numerics.solver.preconditioner.partitioning = partitioning;
    fprintf(1,' done\n');
    simulation = simulation.numerics.prepareW(simulation);
    simulation.numerics.solver.preconditioner.partitioningIdcs = simulation.numerics.partitioningIdcs;
end

%compute all the Mie coefficients
simulation = simulation.computeParallelMieCoefficients;
%Compute scattered field coefficients b_i,n
simulation = simulation.numerics.solver.preconditioner.prepareM(simulation);
%initial field
simulation = simulation.computeInitialFieldCoefficients();
%compute scattered fields
simulation = simulation.computeScatteredFieldCoefficients();

%plot all fields
[x,y] = meshgrid(-3000:50:3000,-3000:50:3000); z=x-x+2000;
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

