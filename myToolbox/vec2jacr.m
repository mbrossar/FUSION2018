function [ J ] = vec2jacr( vec )
tolerance = 1e-12;

if size(vec,1) == 3
    
    phi = vec;
    
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        J = vec2jacrSeries(phi,10);
    else
        axis = phi/ph;

        cph = (1 - cos(ph))/ph;
        sph = sin(ph)/ph;

        J = sph * eye(3) + (1 - sph) * axis * axis' - cph * hat(axis);
    end       
    
elseif size(vec,1) == 6
        
    phi = vec(1:3);
    rho = vec(4:6);
    
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        J = vec2jacrSeries(phi,10);
    else
        Jsmall = vec2jacl( phi );
        Q = vec2Qr( vec );
        J = [ Jsmall Q; zeros(3) Jsmall ];
    end
end
end




