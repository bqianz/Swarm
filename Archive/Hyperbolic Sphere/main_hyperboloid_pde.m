
close all; clear;

%[x,y,z]=sphere(N);

N=12; % N must be even
C = 1/4*pi;

% grid i (N+1) x (N+1), N x N patches
th_max = 1.5;
th_vec = acosh(linspace(1,cosh(th_max),N+1));
phi_vec = (0:N)*2*pi/N;
[phi,th] = meshgrid(phi_vec,th_vec);
[x,y,z] = hyp2cart_unit(th,phi);

% area of each patch
patch_area = 2*pi*(cosh(th_max)-1)/N^2;

up_tol = floor(realmax/N^2/max(1,patch_area));



% grid j N x N
th2_vec = (th_vec(1:end-1) + th_vec(2:end))/2;
phi2_vec = (phi_vec(1:end-1) + phi_vec(2:end))/2;
[phi2,th2] = meshgrid(phi2_vec,th2_vec);
coshth2=cosh(th2);
sinhth2=sinh(th2);

[x2,y2,z2] = hyp2cart_unit(th2,phi2);

%% random initial state
S=size(x);

rho = rand(S)+0.1;


% patch
% rho(3:8,8:12)=0.1;

rho=recalibrate_rho_hyp(rho);
rho(end-2:end,:) = 10^-2;

%% calculate total mass
rho2=rho_jgrid_hyp(rho,th_vec,phi_vec,th2_vec,phi2_vec);
mass = sum(sum(rho2))*patch_area;
final_density=mass/(2*pi);

%% scale rho,rho2
rho = rho/mass; % normalize mass to 1
rho2=rho_jgrid_hyp(rho,th_vec,phi_vec,th2_vec,phi2_vec);
mass = sum(sum(rho2))*patch_area;
final_density=1/(2*pi);

rho(1,:) = final_density;
%% make poles nice
% rho([1:2,end-1:end],:)=final_density;



%% plot initial
hyp = surfl(x,y,z);
hyp.CData=rho;
% caxis([0, max(max(rho))*1.2]);
caxis([0, final_density*2]);
% set(sph, 'FaceAlpha', 1)
shading interp;
% axis equal;
az = 45;
el = 20;
% el = -90;
view(az, el);
colorbar;
% scatter3(x2(:),y2(:),z2(:),'.')
title('t = 0');


%% rho values at j-grid

% fnplt(rho_interp);
% xlabel('th'), ylabel('phi');

% check if interpolation is correct
% rho2=rho_jgrid(rho,th_vec,phi_vec);
% figure(2)
% hyp2 = surfl(x2,y2,z2);
% hyp2.CData=rho2;
% caxis([0, max(max(rho))*1.2]);
% % caxis([0, final_density*2]);
% % set(sph, 'FaceAlpha', 1)
% % shading interp;
% axis square;
% az = 45;
% el = 20;
% % el = -90;
% view(az, el);
% colorbar;
% rho_interp_check=csape({th2_vec,phi2_vec},rho_2);
% rho_check=fnval(rho_interp_check,{th_vec,phi_vec});
% rho_check=recalibrate_rho(rho_check);
% diff=abs(rho-rho_check);
% norm(diff(2:end-1,2:end-1))
% figure(3)
% sph3=surfl(x,y,z);
% sph3.CData=rho_check;
% shading interp;
% axis square;
% az = 45;
% el = 20;
% % el = -90;
% view(az, el);
% colorbar;




%% time-stepping set-up
tmax = 1000;
dt = 10^-3;

drho_dt=zeros(S);

% i-grid
len1=length(th_vec);
len2=length(phi_vec);


pause;
%% begin time stepping

for t=dt:dt:tmax
    
    % space derivative of rho
    [drho_dth, drho_dphi]=drho_dspace_hyp(rho,th_vec,phi_vec,N,S);
    
    
    
    for l=1:len2-1 % copy phi=2*pi from phi=0 later
        phi_i = phi_vec(l);
        
        for k=2:len1
            th_i = th_vec(k);
            sinhth = sinh(th_i);
            coshth = cosh(th_i);
            
            z=coshth*coshth2 - sinhth*sinhth2.*cos(phi_i-phi2);
            dist = acosh(z);
%             
            dzdth = sinhth*coshth2 - coshth*sinhth2*cos(phi_i-phi2);
            
            
            dzdphi = sinhth*sinhth2.*sin(phi_i-phi2);
            ind = dzdphi > 100;
            dzdphi(ind) = 100;
%             
%             dKdz =  -1/(4*pi)./(sinh(dist/2).*cosh(dist/2).*sqrt(z.^2-1) )...
%                 + C*coth(dist) ./sqrt(z.^2-1);

            dKddist = -1./(2*pi*sinh(dist)) + tanh(dist)/(4*pi);
            ddistdth = dzdth./sqrt(z.^2-1);
            ddistdphi = dzdphi./sqrt(z.^-1);
            
           
            mat1 = dKddist.*ddistdth.*rho2;%.*sinth2;
            mat2 =  dKddist.*ddistdphi./sinhth  .*  rho2;%  .*  sinth2;
            
            %             I = mat1>N^2;
            %             mat1(I) = 0;
            %             J = mat2>N^2;
            %             mat2(J)=0;
            
%             I1 = mat1>up_tol;
%             I2 = mat1<-up_tol;
%             mat1(I1) = up_tol;
%             mat1(I2) = -up_tol;
%             
%             J1 = mat2>up_tol;
%             J2 = mat2<-up_tol;
%             mat2(J1) = up_tol;
%             mat2(J2) = -up_tol;
            
            
            temp = rho(k,l)*(final_density - rho(k,l)) + ...
                drho_dth(k,l)*patch_area * sum(sum(  mat1 ));
            
            if k==1 % derivative has no dependence on phi when th = 0, pi
                drho_dt(k,l) = temp;
            else
                drho_dt(k,l) = temp + drho_dphi(k,l)*patch_area./(sinhth^2*...
                    sum(sum(  mat2  )) );
            end
        end
    end
    
    
    
    % make periodic in phi
    drho_dt(:,end) = drho_dt(:,1);
    
    drho_dt(1,:) = 0;
    
    ind = abs(drho_dt) > 10^3;
    drho_dt(ind) = 10^3;
    
    rho = rho + dt*drho_dt;
    
    rho = recalibrate_rho_hyp(rho);
    hyp.CData=rho;
    title(['t = ' num2str(t)]);
    % caxis([0, max(max(rho))*1.2]);
    drawnow;
    
            % rho at j-grid
    rho2=rho_jgrid_hyp(rho,th_vec,phi_vec,th2_vec,phi2_vec);
    
    current_mass = sum(sum(rho2))*patch_area
    
    if norm([drho_dth(:);drho_dphi(:)],inf) < 10^-2/2
        break
    end
    
    
    
end








