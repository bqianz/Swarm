function [drho_dth, drho_dphi]=drho_dspace(rho,th_vec,phi_vec,N,S)
drho_dth = zeros(S);
drho_dphi = zeros(S);

for i=1:N+1
    % for dth: iterate through each column of rho, one column one phi value
    
    p=csape(th_vec,[0,rho(:,i)',0],[1,1]); % neumann boundary
%     mesh=linspace(0,pi,100);
%     interpolated= fnval(p,mesh);
%     plot(th_vec,rho(:,i),'*',mesh,interpolated,'--')
    p_der=fnder(p,1);
    values=ppval(p_der,th_vec);
    values([1,end])=0;
    drho_dth(:,i)=values;
    
    % for dphi: iterate through each row of rho
    p=csape(phi_vec,rho(i,:),[0,0]); % periodic boundary
    p_der=fnder(p,1);
    values=ppval(p_der,phi_vec);
    avg=(values(1)+values(end))/2;
    values([1,end])=avg;
    drho_dphi(i,:) = values;
end
drho_dphi(1,:) = 0;

end