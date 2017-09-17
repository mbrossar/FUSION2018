function [ invJ ] = vec2jaclInv( vec )
tolerance = 1e-12;

if size(vec,1) == 3
    
    phi = vec;
    
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        invJ = vec2jacInvSeries(phi,10);
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        invJ =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   
    
elseif size(vec,1) == 6

    phi = vec(1:3);
    rho = vec(4:6)
    
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        invJ = vec2jaclInvSeries(phi,10);
    else
        invJsmall = vec2jaclInv( phi );
        Q = vec2Q( vec );
        invJ = [ invJsmall -invJsmall*Q*invJsmall; zeros(3) invJsmall ];
    end  
end    
end

