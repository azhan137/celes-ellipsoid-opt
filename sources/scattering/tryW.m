II = eye(120);

newWt = zeros(120);

parfor i = 1:120
    Wi = gather(coupling_matrix_multiply_T(simulation,II(:,i)));
    newWt(:,i) = Wi;
end

figure
imagesc(abs(newWt))
figure
imagesc(abs(Wt))

