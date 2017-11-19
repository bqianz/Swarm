function density = colour_density_hyp(th,phi,mean_short)
num_sets = 100;
N = length(th);
densities = zeros(num_sets,N);
max_th = max(th);

for n = 1:N
    
    th_n = th(n);
    phi_n = phi(n);

    for k = 1:num_sets % trials to average over
        r = (3+rand(1)*2)*mean_short;
        if false % max_th - th_n >= 5*mean_short
            angles = fsolve(@angle_solvers,[0,2,cosh(r),cosh(th_n)]);
            area = 2*angles(1)*2 + (2*pi - 2*angles(3))*cosh(r) + 2*angles(2) - 2*angles(4)-2*pi;
        else
            area = 2*pi*(cosh(r)-1);
        end
        I = acosh(cosh(th_n)*cosh(th) - sinh(th_n)*sinh(th).*cos(phi_n - phi)) < r;
        densities(k,n) = sum(I)/N/area;
    end
end
density = mean(densities);
%density = densities;
end