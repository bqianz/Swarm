function density = colour_density_sph(th,phi)
% generate a set of grids to calculate density and take average of
num_sets = 100;
N = length(th);

densities = zeros(num_sets,N);

grid_scale = pi/12;
for n = 1:N
    
    th_n = th(n);
    phi_n = phi(n);

    for k = 1:num_sets % number of patches to average over
        th_range = max(rand(1)*grid_scale*2,grid_scale);
        phi_range = max(rand(1)*grid_scale*2,grid_scale);
        
        th1 = th_n - (rand(1)*0.8+0.1)*th_range;
        th2 = th1 + th_range;
        
        phi1 = phi_n - (rand(1)*0.8+0.1)*phi_range;
        phi2 = phi1 + phi_range;
        
        phi1 = mod(phi1,2*pi);
        phi2 = mod(phi2,2*pi);
        
        phi3 = mod(phi1+pi,2*pi);
        phi4 = mod(phi2+pi,2*pi);
        
        
        if th1 <= 0
            area = phi_range * (abs( cos(th1) - 1) + abs(1 - cos(th2)));
            
            I = (th==0) + (th>0).*(  (th<=th2).*(phi>=phi1).*(phi<=phi2)...
                + (th<=-th1).*(phi>=phi3).*(phi<=phi4));
            
        elseif th2 >= pi
            area = phi_range * (abs(cos(th1) + 1 ) + abs(-1 - cos(th2)));
            
            I = (th==pi) + (th<pi).*(  (th>=th1).*(phi>=phi1).*(phi<=phi2)...
                + (th>=(2*pi - th2)).*(phi>=phi3).*(phi<=phi4));
            
        else % theta between 0 and pi
            area = phi_range * abs( cos(th1) - cos(th2) );
            I = (th>=th1).*(th<=th2).*(phi>=phi1).*(phi<=phi2);
            
        end
        densities(k,n) = sum(I)/N/area;
    end
end
density = mean(densities);
%density = densities;
end