function [chi,omega_b,a_b,P] = iekfUpdate(chi,omega_b,a_b,...
    P,y,R,ParamFilter,yAmers)
Pi = ParamFilter.Pi;
chiC = ParamFilter.chiC;
RotC = chiC(1:3,1:3);
xC = chiC(1:3,4);
q = length(P);
Rot = chi(1:3,1:3);
x = chi(1:3,5);
NbAmers = length(yAmers);
k = length(y);
R = kron(eye(k/2),R);
H = zeros(k,q);
PosAmers = chi(1:3,6:end);
posAmers = PosAmers(:,yAmers);
z = Pi*( (Rot*RotC)'*(posAmers-kron(x,ones(1,NbAmers))) ...
    - kron(xC,ones(1,NbAmers)));
y_bar = z(1:2,:)./z(3,:);
y_bar = y_bar(:);

for j = 1:length(yAmers)
    z_j = z(:,j); 
    H_j = zeros(3,length(P));
    H_j(:,7:9) = -(Rot*RotC)';
    H_j(:,13+yAmers(j)*3:15+yAmers(j)*3) = (Rot*RotC)';
    delta_h = [1/z_j(3) 0 -z_j(1)/z_j(3)^2;
        0 1/z_j(3) -z_j(2)/z_j(3)^2]*Pi;
    H_j = delta_h*H_j;
    H(2*j-1:2*j,:) = H_j;
end

S = H*P*H' + R;
K = P*H'/S;
P = (eye(q)-K*H)*P;

xibar = K*(y-y_bar); % innovation
omega_b = omega_b + xibar(10:12);
a_b = a_b + xibar(13:15);
xibar = xibar([1:9 16:q]);

% Update mean state
chi = exp_multiSE3(xibar)*chi;
end
