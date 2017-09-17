function [Rot,v,x,PosAmers] = chi2state(chi)
Rot = chi(1:3,1:3);
v = chi(1:3,4);
x = chi(1:3,5);
if size(chi,2) > 5
    PosAmers = chi(1:3,6:end);
else
    PosAmers = [];
end
end
