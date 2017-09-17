function res = logSO3(R)
theta = acos((trace(R)-1)/2) ;
if theta == 0
    res = zeros(3);
else
    res = theta/(2*sin(theta))*(R-R');
end
res = [-res(2,3);res(1,3);-res(1,2)];
end
