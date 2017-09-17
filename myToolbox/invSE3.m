function [ichi] = invSE3(chi)
ichi = eye(size(chi));
R = chi(1:3,1:3);
iR = R';
ichi(1:3,1:3) = iR;
ichi(1:3,4:end) = -iR*chi(1:3,4:end);
end

