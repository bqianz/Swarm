N= 99;
th_i = 0;
phi_i = 0;

th_j = rand(1,N)' * pi;
phi_j = rand(1,N)' * 2*pi;

value = -sin(th_i) * cos(th_j) + cos(th_i) * sin(th_j) * cos(phi_i-phi_j);



