function drho=recalibrate_drho(drho)
% for space derivative of rho only
drho(1,:) = 0;
drho(end,:) = 0;
drho(:,end)=drho(:,1);
end