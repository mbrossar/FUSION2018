function [ invJ ] = vec2jacrInv( vec )
tolerance = 1e-12;

if size(vec,1) == 3
    
    phi = vec;
    
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        invJ = vec2jacrInvSeries(phi,10);
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        invJ =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   
    
elseif size(vec,1) == 6

    rho = vec(1:3);
    phi = vec(4:6);
    
    ph = norm(phi);
    if ph < tolerance;
        % If the angle is small, fall back on the series representation
        invJ = vec2jacInvSeries(phi,10);
    else
        invJsmall = vec2jacrInv( phi );
        Q = vec2Qr( vec );
        invJ = [ invJsmall -invJsmall*Q*invJsmall; zeros(3) invJsmall ];
    end
end    
end

