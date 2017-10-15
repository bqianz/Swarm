function colour_density(th,phi)
global N th_fine phi_fine hyp th_coarse phi_coarse Area n N_fine;
den_fine = zeros(N_fine);
for i = 1:n
    th1 = th_coarse(i);
    th2 = th_coarse(i+1);
    for j = 1:n
        ph1 = phi_coarse(j);
        ph2 = phi_coarse(j+1);
        index1 = and(th>=th1,th<=th2);
        index2 = and(phi>=ph1, phi<=ph2);
        density = sum(index1.*index2) / Area(i,j);
        
        index3 = and(th_fine>=th1,th_fine<=th2);
        index4 = and(phi_fine>=ph1,phi_fine<=ph2);
        den_fine(index3',index4) = density/N;
    end
end
% den = imresize(den,goal_size);
hyp.CData=den_fine;
end