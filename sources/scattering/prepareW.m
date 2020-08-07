function W = prepareW(simul)

M = prepareM(simul);
T = simul.tables.mieCoefficients(simul.input.particles.radiusArrayIndex,:);
invT = inv(diag(T(:)));
W = invT*(eye(simul.tables.nmax*simul.tables.pmax)-M);