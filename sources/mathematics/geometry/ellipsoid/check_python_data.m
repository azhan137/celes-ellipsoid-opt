r_0 = 200;
h_0 = 460;
lambda = 1000;
ns = 2;
ni = 1;

step_size = 0.01;
Ntheta = 100;
lmax = 4;

[J11, J12, J21, J22, dJ11, dJ12, dJ21, dJ22] = compute_J_cyl_celes_posm(lmax,Ntheta,200,460,ni,ns,lambda,3);
[J11j, J12j, J21j, J22j, dJ11j, dJ12j, dJ21j, dJ22j] = compute_J_cyl_celes_posm(lmax,Ntheta,200,460,ni,ns,lambda,1);
