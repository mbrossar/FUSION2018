function [calJl] = xi2calJl(xi)
%Compute Left Jacobian
taille_xi = length(xi);
tolerance = 1e-10;
NbAmers = taille_xi/3-3;
phi = xi(1:3);

ph = norm(phi);
if ph < tolerance
    % If the angle is small, fall back on the series representation
    J = vec2jaclSeries(phi,10);
else
    axis = phi/ph;
    
    cph = (1 - cos(ph))/ph;
    sph = sin(ph)/ph;
    
    J = sph * eye(3) + (1 - sph) * axis * axis' + cph * hat(axis);
    for i = 1:2
        Q = vec2Ql([phi;xi(1+3*i:3+3*i)]);
        calJl(1:3,1+3*i:3+3*i) = Q;
    end
    for i = 1:NbAmers
        Q = vec2Ql([phi;xi(7+3*i:9+3*i)]);
        calJl(1:3,13+3*i:15+3*i) = Q;
    end
end

calJl = kron(eye(NbAmers+5),J);
calJl(10:15,10:15) = eye(6);
end

