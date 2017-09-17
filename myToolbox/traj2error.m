function [error] = traj2error(trajFilter,trajReal)
NbMax = length(trajReal.psi);
error.R = zeros(3,NbMax);
for i = 1:NbMax
    RotFilter = eul2rotm([trajFilter.psi(i+2),trajFilter.theta(i+2),trajFilter.phi(i+2)]);
    RotReal = eul2rotm([trajReal.psi(i+2),trajReal.theta(i+2),trajReal.phi(i+2)]);
    error.R(:,i) = logSO3(RotFilter'*RotReal,'R3');
end
error.v = trajFilter.v(:,1:NbMax)-trajReal.v(:,1:NbMax);
error.x = trajFilter.x(:,1:NbMax)-trajReal.x(:,1:NbMax);
error.omega_b = trajFilter.omega_b(:,1:NbMax)-trajReal.omega_b(:,1:NbMax);
error.a_b = trajFilter.a_b(:,1:NbMax)-trajReal.a_b(:,1:NbMax);

%rms
R = zeros(NbMax,1);
v = zeros(NbMax,1);
x = zeros(NbMax,1);
omega_b = zeros(NbMax,1);
a_b = zeros(NbMax,1);
for i = 1:NbMax
    R(i) = norm(error.R(:,i))/3*180/pi;
    v(i) = norm(error.v(:,i))/3;
    x(i) = norm(error.x(:,i))/3;
    omega_b(i) = norm(error.omega_b(:,i))/3;
    a_b(i) = norm(error.a_b(:,i))/3;
end
error.rms.R = R;
error.rms.v = v;
error.rms.x = x;
error.rms.omega_b = omega_b;
error.rms.a_b = a_b;

end

