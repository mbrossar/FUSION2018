function [calJr] = xi2calJr(xi)
%Compute Right Jacobian
taille_xi = length(xi);
tolerance = 1e-12;
NbAmers = taille_xi/3-3;
phi = xi(1:3);

ph = norm(phi);
if ph < tolerance
    % If the angle is small, fall back on the series representation
    J = vec2jacrSeries(phi,10);
else
    axis = phi/ph;
    
    cph = (1 - cos(ph))/ph;
    sph = sin(ph)/ph;
    
    J = sph * eye(3) + (1 - sph) * axis * axis' - cph * hat(axis);
    for i = 1:2
        Q = vec2Qr([phi;xi(1+3*i:3+3*i)]);
        calJr(1:3,1+3*i:3+3*i) = Q;
    end
    for i = 1:NbAmers
        Q = vec2Qr([phi;xi(7+3*i:9+3*i)]);
        calJr(1:3,13+3*i:15+3*i) = Q;
    end
end

calJr = kron(eye(NbAmers+5),J);
calJr(10:15,10:15) = eye(6);
end
