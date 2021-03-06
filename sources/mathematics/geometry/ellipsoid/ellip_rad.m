





%calculate surface for theta and phi values given a,b, the axes of the
%ellipsoid.

function radius = ellip_rad(theta,phi,a,b,c)

radius = (((cos(phi)/a).^2+(sin(phi)/b).^2).*sin(theta).^2+(cos(theta)/c).^2);
radius = 1./sqrt(radius);
