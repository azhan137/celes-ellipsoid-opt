function imageRightHandSide = compute_image_rhs(simulation,image,E,image_pts)

diff_array = -(image(:)-sum(abs(E).^2,2)).*conj(E);
% diff_array = -abs(image(:)-sum(abs(E).^2,2)).*conj(E);

imageRightHandSide = compute_adjoint_rhs(simulation,diff_array,image_pts);
