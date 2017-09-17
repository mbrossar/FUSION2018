function [xi] = logSE3(chi)
C = chi(1:3,1:3);
r = chi(1:3,4);
[vecs,eigs] = eig(C);
[~,idx] = min(abs(diag(eigs)-1));
a = vecs(:,idx);

phi = acos(1/2*(trace(C)-1));
if phi == 0
    a = zeros(3,1);
    iJ = eye(3);
else
    C1 = cos(phi)*eye(3) + (1-cos(phi))*a*transpose(a) + sin(phi)*[[0,-a(3),a(2)];[a(3),0,-a(1)];[-a(2),a(1),0]];
    C2 = cos(-phi)*eye(3) + (1-cos(-phi))*a*transpose(a) + sin(-phi)*[[0,-a(3),a(2)];[a(3),0,-a(1)];[-a(2),a(1),0]];
    if norm(C2-C)<norm(C1-C)
        phi = -phi;
    end
    iJ = phi/2*cot(phi/2)*eye(3) + (1-phi/2*cot(phi/2))*a*a' - phi/2*[[0,-a(3),a(2)];[a(3),0,-a(1)];[-a(2),a(1),0]];
end
rho = iJ*r;
xi = real([phi*a;rho]);
end
