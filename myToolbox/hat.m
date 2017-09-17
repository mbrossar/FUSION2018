function [ vechat ] = hat( vec )
if size(vec,1) == 3
    vechat = [  0,     -vec(3),  vec(2);
            vec(3),   0    , -vec(1);
           -vec(2),  vec(1),   0    ];    
elseif size(vec,1) == 6
    vechat = [ hat( vec(4:6,1) ) vec(1:3,1); zeros(1,4) ];    
end
end

