function [ J ] = vec2jacrSeries( vec, N )
if size(vec,1) == 3 
    
    J = eye(3);
    pxn = eye(3);
    px = -hat(vec);
    for n = 1:N
        pxn = pxn*px/(n + 1);    
        J = J + pxn;
    end
    
elseif size(vec,1) == 6    
    
    J = eye(6);
    pxn = eye(6);
    px = -curlyhat(vec);
    for n = 1:N
        pxn = pxn*px /(n + 1);    
        J = J + pxn;
    end      
end
end

