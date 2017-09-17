function [points3d, covariance] = myEsimatePosAmers(pointTracks, ...
    camPoses, cameraParams,ParamFilter)
chiC = ParamFilter.chiC;
RotC = chiC(1:3,1:3);
xC = chiC(1:3,4);
R = eye(2);
k = 2*length(camPoses);
q = 6*length(camPoses);
N_aug = q+k;
Rc = chol(kron(eye(k/2),R));

S = zeros(q);
chi = cell(k/2,1);
for i = 1:k/2
    S(6*i-5:6*i,6*i-5:6*i) = camPoses(i).S([1:3 7:9],[1:3 7:9]);
    chi{i} = [camPoses(i).Orientation camPoses(i).Location;0 0 0 1];
end

S_aug = blkdiag(S,Rc);

% scaled unsented transform
W0 = 1-N_aug/3;
Wj = (1-W0)/(2*N_aug);
gamma = sqrt(N_aug/(1-W0));
alpha = 1;
beta = 2;

% Compute transformed measurement
X = zeros(N_aug,2*N_aug+1); %sigma points
Y = zeros(3,2*N_aug+1);
ybar = zeros(3,1); % Measurement mean
for j = 1:2*N_aug+1
    if j == 1
    elseif j < N_aug+2
        X(:,j) = gamma*S_aug(j-1,:)';
    else
        X(:,j) = -gamma*S_aug(j-N_aug-1,:)';
    end
    camPoses_j = camPoses;
    pointTracks_j = pointTracks;
    xi_j = X(1:q,j);
    v_j = X(q+1:N_aug,j);
    for i = 1:k/2
        chi_j = expSE3(xi_j(6*i-5:6*i))*chi{i};
        camPoses_j(i).Orientation = RotC'*chi_j(1:3,1:3)';
        camPoses_j(i).Location = chi_j(1:3,4) +chi_j(1:3,1:3)*RotC*xC;
        pointTracks_j.Points(i,:) = pointTracks.Points(i,:) + v_j(2*i-1:2*i)';
    end
    
    Y(:,j) = EsimatePosAmers(pointTracks_j,camPoses_j,cameraParams)';
    if j == 1
        ybar = ybar + W0*Y(:,j);
    else
        ybar = ybar + Wj*Y(:,j);
    end
end

Y(:,1) = sqrt(abs(W0+(1-alpha^2+beta)))*(Y(:,1)-ybar);
YY = sqrt(Wj)*(Y(:,2:2*N_aug+1)-ybar*ones(1,2*N_aug));
[~,Rs] = qr(YY');
Ss = Rs(1:3,1:3);
[Sy,~] = cholupdate(Ss,Y(:,1),'-'); % Sy'*Sy = Pyy

Pxy = zeros(q,3);
for j = 2:2*N_aug+1
    Pxy = Pxy + Wj*X(1:q,j)*(Y(:,j)-ybar)';
end

points3d = ybar;
S = S(end-5:end,end-5:end);
covariance = [S'*S Pxy(end-5:end,:); Pxy(end-5:end,:)' Sy'*Sy];
end