%%calculate derivative of specified associated legendre polynomial. 
%d/dtheta[p_{l^m}(ct)]=[l*ct*p_{l^m}(ct)-(l+m)*p_{l-1^m}(ct)]/sqrt(1-ct^2)
%plm: plm table in cell
%l: orbital index non-negative integer
%m: azimuthal index integer
%ct: argument of legendre polynomial

%returns derivative

function deriv_legendre = assoc_legendre_deriv(plm,l,m,ct)

%p_{l^-m} = (-1)^m*(l-m)!/(l+m)!*p_{l^m}
pos_to_neg = (-1).^m*factorial(l-m)/factorial(l+m);

%p_{l^m} = 0 for m > l
if (abs(m) <= l-1)
    deriv_legendre = (l*ct.*plm{l+1,abs(m)+1}-(l+m)*plm{l,abs(m)+1})./sqrt(1-ct.^2);
else
    deriv_legendre = (l*ct.*plm{l+1,abs(m)+1})./sqrt(1-ct.^2);
end
%convert m -> -m if necessary
if (m < 0)
    deriv_legendre = pos_to_neg*deriv_legendre;
end

% old code
% if (m >= 0 && abs(m) <= l-1)
%     deriv_legendre = (l*ct.*plm{l+1,m+1}-(l+m)*plm{l,m+1})./sqrt(1-ct.^2);
% elseif (m >= 0 && abs(m) > l-1)
%     deriv_legendre = (l*ct.*plm{l+1,m+1})./sqrt(1-ct.^2);
% elseif (m < 0 && abs(m) <= l-1)
%     deriv_legendre = pos_to_neg*(l*ct.*plm{l+1,m+1}-(l+m)*plm{l,m+1})./sqrt(1-ct.^2);
% elseif (m < 0 && abs(m) > l-1)
%     deriv_legendre = pos_to_neg*(l*ct.*plm{l+1,m+1})./sqrt(1-ct.^2);
% else
%     disp(l);
%     fprintf('.\n');
%     disp(m);
% end


