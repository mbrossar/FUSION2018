function [ invJ ] = vec2jacrInvSeries( vec, N )

if size(vec,1) == 3 
    
    invJ = eye(3);
    pxn = eye(3);
    px = -hat(vec);
    for n = 1:N
        pxn = pxn * px/n;
        invJ = invJ + bernoullinumber(n) * pxn;
    end
    
elseif size(vec,1) == 6    
    
    invJ = eye(6);
    pxn = eye(6);
    px = -curlyhat(vec);
    for n = 1:N
        pxn = pxn * px/n;
        invJ = invJ + bernoullinumber(n) * pxn;
    end
end
end

