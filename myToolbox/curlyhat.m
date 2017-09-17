function [veccurlyhat] = curlyhat(vec)
phihat = hat( vec(4:6) );
veccurlyhat = [ phihat hat(vec(1:3)); zeros(3) phihat ];
end
