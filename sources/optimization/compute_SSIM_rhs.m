%%derivative of the SSIM

function imageRightHandSide = compute_SSIM_rhs(simulation,img_ref,mean_ref,var_ref,E_iter,mean_iter,var_iter,covar,img_pts)
c1 = 0.0001;
c2 = 0.0009;

denom = (mean_ref^2+mean_iter^2+c1)*(var_ref+var_iter+c2);
N = length(img_pts(:,1));
dmean = ones(N,1)/N;
dvar = 2/N*(sum(abs(E_iter).^2,2)-mean_iter).*(1-dmean);
dcovar = 1/N*(img_ref-mean_ref).*(1-dmean);

num1 = 2*mean_ref*dmean*(2*covar+c2)+(2*mean_iter*mean_ref+c1)*2*dcovar;
num2 = (2*mean_ref*mean_iter+c1)*(2*covar+c2)*(2*mean_ref*dmean*(var_ref+var_iter+c1)+2*sqrt(var_iter)*dvar*(mean_iter^2+mean_ref^2+c1));

dSSIM = num1/denom-num2/denom^2;

diff_array = dSSIM.*conj(E_iter);

imageRightHandSide = compute_adjoint_rhs(simulation,diff_array,img_pts);