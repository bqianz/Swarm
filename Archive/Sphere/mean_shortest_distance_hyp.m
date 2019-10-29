function [mean_nearest, nearest_dist,farthest_dist] = mean_shortest_distance_hyp(th,phi)
N = length(th(:));
nearest_dist=zeros(1,N);
farthest_dist=zeros(1,N);

for i = 1:N
    th1 = th(i);
    phi1 = phi(i);
    temp = sort(acosh(cosh(th1)*cosh(th) - sinh(th1)*sinh(th).*cos(phi1 - phi)));
    nearest_dist(i) = temp(2);
    farthest_dist(i) = temp(end);
end
mean_nearest = mean(nearest_dist);
end
