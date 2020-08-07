function imageRightHandSide = compute_contrast_max_rhs(simulation,E_max,image_pts)
% tis function computes the term
% dFOM/dE = dFOM/dE_max+dFOM/dE_min for a single polarization
% the figure of merit given is 
% FOM = [(I_TE_MAX-I_TE_MIN)+(I_TM_MAX-I_TM_MIN)]
diff_array = conj(E_max);
imageRightHandSide = compute_adjoint_rhs(simulation,diff_array,image_pts);