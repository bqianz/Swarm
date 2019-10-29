function [x,y,z] = hyp2cart_unit(th,phi)

sinhth=sinh(th);
x = sinhth .* cos(phi);
y = sinhth .* sin(phi);
z = cosh(th);
