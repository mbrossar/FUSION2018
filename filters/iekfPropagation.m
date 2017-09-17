function [Rot,v,x,PosAmers,P] = iekfPropagation(dt,Rot,v,x,omega_b,a_b,...
    PosAmers,P,omega,acc,Q,g)
NbAmers = size(PosAmers,2);

% state propagation
Rot = Rot*expSO3((omega-omega_b)*dt);
v = v+(Rot*(acc-a_b)+g)*dt;
x = x+v*dt;

% covariance propagation
F = eye(size(P));
F(4:6,1:3) = vecto(g)*dt;
F(7:9,1:3) = vecto(g)*dt*dt;
F(7:9,4:6) = eye(3)*dt;
F(1:3,10:12) = -Rot*dt;
F(4:6,10:12) = -vecto(v)*Rot*dt;
F(7:9,10:12) = -vecto(x)*Rot*dt;
for i = 1:NbAmers
    posAmers_i = PosAmers(:,i);
    F(13+3*i:15+3*i,10:12) = -vecto(posAmers_i)*Rot*dt;
end
F(4:6,13:15) = -Rot*dt;
F(7:9,13:15) = -Rot*dt*dt;

G = zeros(size(P,1),size(Q,1));
G(1:3,1:3) = Rot;
G(4:6,1:3) = vecto(v)*Rot;
G(7:9,1:3) = vecto(x)*Rot;
G(4:6,4:6) = Rot;
G(7:9,4:6) = Rot*dt;
G(10:15,7:12) = eye(6);
G(1:3,7:9) = Rot*dt;
G(4:6,7:9) = vecto(v)*Rot*dt*dt;
G(7:9,7:9) = vecto(x)*Rot*dt*dt*dt;
G(4:6,10:12) = Rot*dt*dt;
G(7:9,10:12) = Rot*dt*dt*dt;

P = F*P*F' + G*(Q*dt)*G'*dt;
end
