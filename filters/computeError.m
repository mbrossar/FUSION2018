function [error] = computeError(trajFilter,trajReal,i)
errorR = zeros(i-2,1);
errorX = zeros(i-2,1);
for j = 2:i-1
    RotFilter = squeeze(trajFilter.Rot(:,:,j));
    RotReal = eul2rotm([trajReal.psi(j),trajReal.theta(j),trajReal.phi(j)]);
    errorR(j-1) = norm(logSO3(RotFilter'*RotReal));
    errorX(j-1) = norm(trajFilter.x(:,j)-trajReal.x(:,j));
end
errorR = 180/pi*errorR/3; %for RMSE
errorX = errorX/3; %for RMSE
error.errorR = errorR;
error.errorX = errorX;
end