%% plot sphere
close all; clear;
R=1;


%[x,y,z]=sphere(N);

N=20; % N must be even

% grid i (N+1) x (N+1), N x N patches
th_vec = acos((N/2:-1:-N/2)*2/N);
phi_vec = (0:N)*2*pi/N;
[phi,th] = meshgrid(phi_vec,th_vec);
[x,y,z] = sph2cart(phi,pi/2-th,1);

% area of each patch
patch_area = 4*pi/N^2;

up_tol = floor(realmax/N^2/max(1,patch_area));


% grid j N x N
th2_vec = (th_vec(1:end-1) + th_vec(2:end))/2;
phi2_vec = (phi_vec(1:end-1) + phi_vec(2:end))/2;
[phi2,th2] = meshgrid(phi2_vec,th2_vec);
costh2=cos(th2);
sinth2=sin(th2);

[x2,y2,z2] = sph2cart(phi2,pi/2-th2,1+eps);

%% random initial state
S=size(x);

rho = rand(S)+0.1;

% patch
% rho(3:8,8:12)=0.1;

rho=recalibrate_rho(rho);

% zero at poles
% rho(1,:)=0;
% rho(end,:)=0;

%% calculate total mass
rho2=rho_jgrid(rho,th_vec,phi_vec,th2_vec,phi2_vec);
mass = sum(sum(rho2))*patch_area;
final_density=mass/(4*pi);

%% scale rho,rho2
rho = rho/mass; % normalize mass to 1
rho2=rho_jgrid(rho,th_vec,phi_vec,th2_vec,phi2_vec);
mass = sum(sum(rho2))*patch_area;
final_density=1/(4*pi);

%% make poles nice
% rho([1:2,end-1:end],:)=final_density;



%% plot initial
sph = surfl(x,y,z);
sph.CData=rho;
caxis([0, final_density*2]);
% set(sph, 'FaceAlpha', 1)
shading interp;
axis square;
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
% sph2 = surfl(x2,y2,z2);
% sph2.CData=rho_2;
% % set(sph, 'FaceAlpha', 1)
% shading interp;
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
dt = 10^-1;

drho_dt=zeros(S);

% i-grid
len1=length(th_vec);
len2=length(phi_vec);


pause;
%% begin time stepping

for t=dt:dt:tmax
    
    % rho at j-grid
    rho2=rho_jgrid(rho,th_vec,phi_vec,th2_vec,phi2_vec);
    
    % space derivative of rho
    [drho_dth, drho_dphi]=drho_dspace(rho,th_vec,phi_vec,N,S);
    
    
    
    for l=1:len2-1 % copy phi=2*pi from phi=0 later
        phi_i = phi_vec(l);
        
        for k=1:len1
            th_i = th_vec(k);
            
            % d is a j-grid matrix corresponding to each i-grid point
            cosdist = cos(th_i)*costh2 + sin(th_i)*sinth2.*cos(phi_i-phi2);
            
            dK_ddist=-cot( acos(cosdist)/2 )./(4*pi);
            
            dacos = -1./sqrt(1-cosdist.^2);
            
            ddist_dth =  ( -sin(th_i)*costh2 + cos(th_i)*sinth2.*cos(phi_i-phi2) ).*dacos;
            
            ddist_dphi = (-sin(th_i)*sinth2.*sin(th_i-th2)).*dacos;
            
%             surf(th2,phi2,dK_ddist.*ddist_dth.*rho2.*sinth2);
%             drawnow;
%             pause
            
            mat1 = dK_ddist.*ddist_dth.*rho2;%.*sinth2;
            mat2 = dK_ddist  .*  ddist_dphi  .*  rho2;%  .*  sinth2;
            
%             I = mat1>N^2;
%             mat1(I) = 0;
%             J = mat2>N^2;
%             mat2(J)=0;

            I1 = mat1>up_tol;
            I2 = mat1<-up_tol;
            mat1(I1) = up_tol;
            mat1(I2) = -up_tol;
            
            J1 = mat2>up_tol;
            J2 = mat2<-up_tol;
            mat2(J1) = up_tol;
            mat2(J2) = -up_tol;
            

            temp = rho(k,l)*(mass/(4*pi) - rho(k,l)) + ...
                drho_dth(k,l)*patch_area * sum(sum(  mat1 ));
            
            if or(k==1,k==len1) % derivative has no dependence on phi when th = 0, pi
                drho_dt(k,l) = temp;
            else
                drho_dt(k,l) = temp + drho_dphi(k,l)*patch_area/(sin(th_i))^2*...
                    sum(sum(  mat2  ));
            end
        end
    end
    
   
    
    temp = (drho_dt(:,end) + drho_dt(:,1))./2;
    drho_dt(:,end) = temp;
    drho_dt(:,1) = temp;
    
    % unify values corresponding to different phi at th=0
    drho_dt([1,end],:) = 0;
    
    rho = rho + dt*drho_dt;
    
    rho = recalibrate_rho(rho);
    sph.CData=rho;
    title(['t = ' num2str(t)]);
    drawnow;
    
%     norm(drho_dth)
%     norm(drho_dphi)

mass = sum(sum(rho2))*patch_area

    if norm([drho_dth(:);drho_dphi(:)],inf) < 10^-2/2
        break
    end
    

    
end








