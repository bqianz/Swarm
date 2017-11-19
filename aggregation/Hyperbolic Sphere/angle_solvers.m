function F = angle_solvers(angles)
% angles = [alpha, beta, gamma, delta]
F(1) = cos(angles(1))./sin(angles(2)) - cos(angles(3))./sin(angles(4));
F(2) = cot(angles(1)).*cot(angles(2));
F(3) = cot(angles(3)).*cot(angles(4));
F(4) = ( cos(angles(2)-angles(4)) + cos(angles(1)).*cos(pi-angles(3)) )...
    ./ ( sin(angles(1)).*sin(pi-angles(3)) );
end