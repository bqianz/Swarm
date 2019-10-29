function rho=recalibrate_rho(rho)
% for rho i-grid, and rho time derivative

% make periodic in phi
rho(:,end) = rho(:,1);

% col = (rho(:,1)+rho(:,end))/2;
% rho(:,1) = col;
% rho(:,end) = col;

% unify values corresponding to different phi at th=0,pi
rho(1,:) = mean(rho(1,:));
rho(end,:) = mean(rho(end,:));
end