function [chi,v,PosAmers,omega_b,a_b,S,xidot] = ukfRefUpdate(chi,v,omega_b,a_b,...
    S,y,param,R,ParamFilter,PosAmers,xidot)
param.Pi = ParamFilter.Pi;
param.chiC = ParamFilter.chiC;

k = length(y);
q = length(S);
N_aug = q+k;
Rc = chol(kron(eye(k/2),R));
S_aug = blkdiag(S,Rc);

% scaled unsented transform
W0 = 1-N_aug/3;
Wj = (1-W0)/(2*N_aug);
gamma = sqrt(N_aug/(1-W0));
alpha = 1;
beta = 2;

% Compute transformed measurement
X = gamma*[zeros(N_aug,1) S_aug' -S_aug'];% sigma-points
Y = zeros(k,2*N_aug+1);
Y(:,1) = h(chi,zeros(q-9,1),param,zeros(N_aug-q,1));
for j = 2:2*N_aug+1
    xi_j = X([1:3 7:9 16:q],j);
    v_j = X(q+1:N_aug,j);
    Y(:,j) = h(chi,xi_j,param,v_j);
end
ybar = W0*Y(:,1) + Wj*sum(Y(:,2:end),2);% Measurement mean
Y(:,1) = sqrt(abs(W0+(1-alpha^2+beta)))*(Y(:,1)-ybar);
YY = sqrt(Wj)*(Y(:,2:2*N_aug+1)-ybar*ones(1,2*N_aug));
[~,Rs] = qr(YY');
Ss = Rs(1:k,1:k);
[Sy,~] = cholupdate(Ss,Y(:,1),'-'); % Sy'*Sy = Pyy

Pxy = zeros(q,k);
for j = 2:2*N_aug+1
    Pxy = Pxy + Wj*X(1:q,j)*(Y(:,j)-ybar)';
end

K = Pxy*Sy^-1*Sy'^-1; % Gain

xibar = K*(y-ybar);
omega_b = omega_b + xibar(10:12);
a_b = a_b + xibar(13:15);
xibar = xibar([1:9 16:q]);

% Covariance update
A = K*Sy';
for n = 1:k
    [S,~] = cholupdate(S,A(:,n),'-');
end
PosAmers = PosAmers + reshape(xibar(10:end),[3 length(xibar(10:end))/3]);
% Update mean state
v = v + xibar(4:6);
chi = [chi(1:3,1:3) chi(1:3,5);0 0 0 1]*expSE3(xibar([1:3 7:9]));

% Parallel transport
B = vecto(xidot(1:3));
C = vecto(xidot(4:6));
alphaB = norm(xidot(1:3));
alphaC = 1/2*norm(xidot(4:6));
B = (eye(3) +1/(alphaC^2)*(1-cos(alphaB))*C + 1/(alphaC^3)*C^2)*B;
C = eye(3) + 1/(alphaC^2)*(1-cos(alphaC))*C^2 + sin(alphaB)/alphaC*C;
expA = [C zeros(3);
    B eye(3)];
P = S'*S;
P([1:3 7:9],[1:3 7:9]) = expA*P([1:3 7:9],[1:3 7:9])*expA';
P([1:3 7:9],[4:6,10:end]) = expA*P([1:3 7:9],[4:6,10:end]);
P([4:6,10:end],[1:3 7:9]) = P([1:3 7:9],[4:6,10:end])';
S = chol(P);
xidot = zeros(6,1);
end

%--------------------------------------------------------------------------
function y = h(chi,xi,param,v)
Pi = param.Pi;
chiC = param.chiC;
RotC = chiC(1:3,1:3);
xC = chiC(1:3,4);

yAmers = param.yAmers;
NbAmers = length(yAmers);

chi_j = chi([1:3 5],[1:3 5])*expSE3(xi(1:6));
Rot = chi_j(1:3,1:3);
x = chi_j(1:3,4);
PosAmers = chi(1:3,6:end) + reshape(xi(7:end),[3 length(xi(7:end))/3]);
posAmers = PosAmers(:,yAmers);
z = Pi*( (Rot*RotC)'*(posAmers-kron(x,ones(1,NbAmers))) ...
    - kron(xC,ones(1,NbAmers)));
y = z(1:2,:)./z(3,:);
y = y(:) + v;
end
