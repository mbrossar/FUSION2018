function [S] = ukfRefPropagation(dt,chi,chiAnt,vAnt,omega_b,a_b,S,omega,acc,Qc,g)
q = length(S); % state size
S_aug = blkdiag(S,Qc);
N_aug = length(S_aug);

%scaled unsented transform
W0 = 1-N_aug/3;
Wj = (1-W0)/(2*N_aug);
gamma = sqrt(N_aug/(1-W0));

ichi = invSE3(chi);
omega = omega-omega_b; % unbiased inputs
acc = acc-a_b;
X = gamma*[zeros(N_aug,1) S_aug' -S_aug'];% sigma-points
X(10:15,:) = X(10:15,:)+X(q+7:q+12,:)*dt;
for j = 2:2*N_aug+1
    xi_j = X([1:3 7:9],j); %do not take bias and landmarks
    xi_v = X(4:6,j);
    w_j = X(q+1:N_aug,j);
    omega_bj = X(10:12,j);
    a_bj = X(13:15,j);
    chi_j = chiAnt*expSE3(xi_j);
    Rot_j = chi_j(1:3,1:3)*expSO3((omega+w_j(1:3)-omega_bj)*dt);
    v_j = vAnt + xi_v + (Rot_j*(acc+w_j(4:6)-a_bj)+g)*dt;
    x_j = chi_j(1:3,4)+v_j*dt;
    chi_j = [Rot_j x_j; zeros(1,3) 1];
    Xi_j = ichi*chi_j; %can be more time efficient
    X([1:3 7:9],j) = logSE3(Xi_j); %propagated sigma points
    X(4:6,j) = v_j - vAnt;
end
X = sqrt(Wj)*X;
[~,Rs] = qr(X(1:q,2:2*N_aug+1)');
S = Rs(1:q,1:q);
end
