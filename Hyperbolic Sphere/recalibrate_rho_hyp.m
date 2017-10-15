function rho=recalibrate_rho_hyp(rho)
% for rho i-grid, and rho time derivative

% make periodic in phi
temp = (rho(:,end) + rho(:,1))./2;
rho(:,end) = temp;
rho(:,1) = temp;

% unify values corresponding to different phi at th=0
rho(1,:) = mean(rho(1,:));

I = (rho<=0);
rho(I) = 10^-2;

rho(end,:) = 10^-2;
end