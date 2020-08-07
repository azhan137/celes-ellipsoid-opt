addpath(genpath('.'))

simulation = celes_simulation;
particles = celes_particles;
input = celes_input;
tables = celes_tables;
initialField = celes_initialField;
numerics = celes_numerics;

lmax = 6;
numParts = 600;

radius = linspace(100,1000,numParts)';
position = rand(numParts,3);
refractiveIndex = ones(1,numParts)*1.52;

particles.positionArray = position;
particles.radiusArray = radius;
particles.refractiveIndexArray = refractiveIndex;

% polar angle of incoming beam/wave, in radians (for Gaussian beams, 
% only 0 and pi are currently possible)
initialField.polarAngle = 0;

% azimuthal angle of incoming beam/wave, in radians
initialField.azimuthalAngle = 0;

% polarization of incoming beam/wave ('TE' or 'TM')
initialField.polarization = 'TE';

% width of beam waist (use 0 or inf for plane wave)
initialField.beamWidth = 0;

% vacuum wavelength (same unit as particle positions and radius)
input.wavelength = 1000;
input.mediumRefractiveIndex = 1;

numerics.lmax = lmax;

input.initialField = initialField;
simulation.numerics = numerics;
simulation.tables = celes_tables;
input.particles = particles;
simulation.input = input;
simulation.tables.pmax = simulation.input.particles.number;
simulation.tables.nmax = simulation.numerics.nmax;

nmax = simulation.tables.nmax;
singleMie = tic;
simulation = simulation.computeParallelMieCoefficients;
timeSingle = toc(singleMie);
% parMie = tic;
% simulation = simulation.computeParallelMieCoefficients;
% timePar = toc(parMie);

mieCoeff = simulation.tables.mieCoefficients;
% figure
% subplot(2,1,1)
% imagesc(radius,1:100,real(mieCoeff))
% title('real mieCoeff');
% colorbar
% subplot(2,1,2)
% imagesc(radius,1:100,imag(mieCoeff))
% title('imag mieCoeff');
% colorbar

% figure
% subplot(2,1,1)
% imagesc(1:nmax,radius,abs(mieCoeff));
% title('abs mieCoeff');
% colorbar
% subplot(2,1,2)
% imagesc(1:nmax,radius,angle(mieCoeff));
% title('phase mieCoeff');
% colorbar

figure
imagesc(1:nmax,radius,abs(mieCoeff));
title('abs mieCoeff');
colormap('gray')
colorbar
