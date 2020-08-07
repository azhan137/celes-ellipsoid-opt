function grad_I = compute_intensity_gradient(simulation,points)

E = compute_scattered_field_opt(simulation,points);
numGradient = simulation.input.particles.number;
grad_I = simulation.numerics.deviceArray(zeros(numGradient,1,'single'));

for g_i = 1:numGradient
    grad_E = compute_grad_scattered_field_opt(simulation,points,g_i);
    grad_I(g_i) = 2*gather(real(sum(grad_E.*conj(E))));
end