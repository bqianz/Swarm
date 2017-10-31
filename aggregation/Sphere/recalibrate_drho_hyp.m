function drho=recalibrate_drho_hyp(drho)
% for space derivative of rho only
drho(1,:) = 0;

drho(:,end)=drho(:,1);
end