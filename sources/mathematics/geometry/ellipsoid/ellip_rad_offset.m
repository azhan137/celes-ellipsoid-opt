%calculate surface for theta and phi values given a,b, the axes of the
%ellipsoid.

function radius = ellip_rad_offset(theta,phi,phi_offset,a,b,c)

radius = (((cos(phi-phi_offset)/a).^2+(sin(phi-phi_offset)/b).^2).*sin(theta).^2+(cos(theta)/c).^2);
radius = 1./sqrt(radius);
