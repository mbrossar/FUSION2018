function chi = expSE3(xi)
theta = norm(xi(1:3));
if(theta == 0)
    chi = [[eye(3),xi(4:6)];[zeros(1,3),1]] ;
else
    Xi = [[vecto(xi(1:3)),xi(4:6)];zeros(1,4)] ;
    chi = eye(4) + Xi + 1/theta^2*(1-cos(theta))*Xi^2 + 1/theta^3*(theta-sin(theta))*Xi^3;
end
end