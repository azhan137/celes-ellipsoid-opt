%
%for computing adjoint derivative
%%FoM = (sum_i(I_0(r_i)-I_a(r_i))^2+sum_j(I_1(r_j)-I_b(r_j))^2))^2
%%adjoint derivative is
%see adjoint image derivative for work
%simulation1: celes_simulation object at wavelength 1
%simulation2: celes_simulation object at wavelength 2
%E1: target electric field at w1 for image
%points1: coordinates of electric field at w1
%E2: target electric field at w2 for image 
%points2: coordinates of electric field at w2


function adjointDerivative = compute_adjoint_chromatic_image_derivative(simulation1,simulation2,E1,E2,points1,points2)
    N1 = compute_bessel_value(simulation1,points1);
    N2 = compute_bessel_value(simulation2,points2);
    
    image_diff1 = image1-sum(E1.*conj(E1));
    image_diff2 = image2-sum(E2.*conj(E2));
    
    adjDeriv1 = image_diff1.*sum(conj(E1).*N1,2);
    adjDeriv2 = image_diff2.*sum(conj(E2).*N2,2);
    
    prefactorSum = sum(image_diff1.^2)+sum(image_diff2.^2);
    
    adjointDerivative = prefactorSum*(squeeze(sum(adjDeriv1,1))+squeeze(sum(adjDeriv2,1)));
end