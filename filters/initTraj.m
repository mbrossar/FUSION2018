function [traj] = initTraj(NbSteps)
traj.Rot = zeros(3,3,NbSteps);
traj.phi = zeros(1,NbSteps);
traj.theta = zeros(1,NbSteps);
traj.psi = zeros(1,NbSteps);
traj.v = zeros(3,NbSteps);
traj.x = zeros(3,NbSteps);
traj.omega_b = zeros(3,NbSteps);
traj.a_b = zeros(3,NbSteps);
end

