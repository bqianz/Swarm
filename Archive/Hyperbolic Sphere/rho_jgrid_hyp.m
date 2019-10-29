function rho2=rho_jgrid_hyp(rho,th_vec,phi_vec,th2_vec,phi2_vec)
len2=length(phi_vec);
rho_clamp_cond = [zeros(1,len2);rho;zeros(1,len2)];
rho_interp = csape({th_vec,phi_vec},rho_clamp_cond,{'clamped','periodic'});
rho2 = fnval(rho_interp,{th2_vec,phi2_vec});

% rho_interp = csape({th_vec,phi_vec},rho);
% rho2 = fnval(rho_interp,{th2_vec,phi2_vec});
end