function density = colour_density_sph(th,phi,mean_short)
num_sets = 100;
N = length(th);
densities = zeros(num_sets,N);

for n = 1:N
    
    th_n = th(n);
    phi_n = phi(n);

    for k = 1:num_sets % trials to average over
        r = (rand(1)+1)*3*mean_short;
        area = 2*pi*(1 - cos(r));
        I = acos(cos(th_n)*cos(th) + sin(th_n)*sin(th).*cos(phi_n - phi)) < r;
        sum(I)
        densities(k,n) = sum(I)/N/area;
    end
end
density = mean(densities);
%density = densities;
end