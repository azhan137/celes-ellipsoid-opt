phi_lims = [0,pi];
Nphi = 30;
[phi_l,w_phi_l] = generate_gauss_weights_abscissae(Nphi,phi_lims(1),phi_lims(2));

m1 = -6:1:6;
m2 = -6:1:6;

numeric_integrals = zeros(length(m1),length(m2));
analytic_integrals = numeric_integrals;

for i = 1:length(m1)
    for j = 1:length(m2)
        phi_exp = exp(1i*(m1(i)-m2(j))*phi_l);
        numeric_integrals(i,j) = sum(w_phi_l.*phi_exp);
        if m1(i) ~= m2(j)
            analytic_integrals(i,j) = -1i*(exp(1i*(m1(i)-m2(j))*pi)-1)/(m1(i)-m2(j));
        else
            analytic_integrals(i,j) = pi;
        end
    end
end

figure
subplot(3,1,1)
imagesc(abs(numeric_integrals))
colorbar
title('numeric')
subplot(3,1,2)
imagesc(abs(analytic_integrals))
colorbar
title('analytic')
subplot(3,1,3)
imagesc(abs(numeric_integrals-analytic_integrals))
colorbar
title('diff')

figure
imagesc(abs(numeric_integrals./analytic_integrals))
colorbar