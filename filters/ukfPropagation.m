function [S] = ukfPropagation(dt,Rot,RotAnt,vAnt,omega_b,a_b,S,omega,acc,Qc,g)
q = length(S);
S_aug = blkdiag(S,Qc);
N_aug = length(S_aug);

% scaled unsented transform
W0 = 1-N_aug/3;
Wj = (1-W0)/(2*N_aug);
gamma = sqrt(N_aug/(1-W0));

ichi = Rot';
omega = omega-omega_b; % unbiased inputs
acc = acc-a_b;
X = gamma*[zeros(N_aug,1) S_aug' -S_aug'];% sigma-points
X(10:15,:) = X(10:15,:)+X(q+7:q+12,:)*dt;
for j = 2:2*N_aug+1
    xi_j = X(1:3,j);
    xi_v = X(4:6,j);
    xi_x = X(7:9,j);
    w_j = X(q+1:N_aug,j);
    omega_bj = X(10:12,j);
    a_bj = X(13:15,j);
    Rot_j = RotAnt*expSO3(xi_j)*expSO3((omega+w_j(1:3)-omega_bj)*dt);
    v_j = vAnt + xi_v + (Rot_j*(acc+w_j(4:6)-a_bj)+g)*dt;
    x_j = xi_x + v_j*dt;
    Xi_j = ichi*Rot_j;
    X(1:3,j) = logSO3(Xi_j); % propagated sigma points
    X(4:6,j) = v_j - vAnt;
    X(7:9,j) = x_j + (v_j-vAnt)*dt;
end
X = sqrt(Wj)*X;
[~,Rs] = qr(X(1:q,2:2*N_aug+1)');
S = Rs(1:q,1:q);
end
