%%
%for computing adjoint
%dFoM/dR_j = dFoM/dE_k*dE_k/db_i,n*db_i,n/dR_j
%adjointDerivative is 
%dFoM/dE_k*dE_k/db_i,n=sum_a(2*(I_0(r_a)-I_j(r_a))conj(E_k(r_a))*N_i,n(r_a))

%simulation is a celes_simulation object
%image is an m x n array of real doubles to be optimized toward
%points is an 3 x m*n array of points corresponding to image

function adjointDerivative = compute_adjoint_image_derivative(simulation,image,points)

conjE = conj(compute_scattered_field_opt(simulation,points));
%conjE = 2*real(conjE);
N = compute_bessel_value(simulation,points);

image_k = sum(conjE.*conj(conjE),2);

adjointDerivative = conjE.*N;
adjointDerivative = squeeze(sum(adjointDerivative,2));

image_a = image-image_k;

adjointDerivative = image_a.*adjointDerivative;
adjointDerivative = squeeze(sum(adjointDerivative,1));
end