function [traj] = updateTraj(traj,Rot,v,x,omega_b,a_b,i)
%record trajectory
angles = rotm2eul(Rot);
traj.Rot(:,:,i) = Rot;
traj.phi(i) = angles(3);
traj.theta(i) = angles(2);
traj.psi(i) = angles(1);
traj.v(:,i) = v;
traj.x(:,i) = x;
traj.omega_b(:,i) = omega_b;
traj.a_b(:,i) = a_b;
end

