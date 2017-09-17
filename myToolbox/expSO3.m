function chi = expSO3(xi)
theta = norm(xi);
if(theta == 0)
    chi = eye(3) ;
else
    Aksi = [[0,-xi(3),xi(2)];[xi(3),0,-xi(1)];[-xi(2),xi(1),0]] ;
    chi = eye(3) + sin(theta)/theta*Aksi + 2*sin(theta/2)^2/theta^2*Aksi^2;
end
end