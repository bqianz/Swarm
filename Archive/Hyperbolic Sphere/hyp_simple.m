%% 2-sphere: potential giving uniform density distribution at equilibrium
clear, close all;
global N th_fine phi_fine hyp th_coarse phi_coarse Area n N_fine;

%% Radius and number of particles


N = 10;
% attractive constant C
C = 1/(4*pi);
q = 1;

%% time-step
dt=10^-2;
tmax=10^8;
step_num=ceil(tmax/dt);
% dt = tmax/step_num;
% frame_skip = floor(step_num/200);
frame_skip = 100;
%% initial condition and plotting it

% regular
% phi=rand(N,1)*2*pi;
% th=rand(N,1)*2;

% tall far patch
% phi=rand(N,1)*pi/10;
% th=rand(N,1) + 3;

% far patch of density similar to steady state density
% phi=rand(N,1)*pi/2;
% th=rand(N,1)*0.3358 + 3;



% two far patches
% N1 = round(N/2);
% N2 = N - N1;
% 
% phi = zeros(N,1);
% th = zeros(N,1);
% 
% phi(1:N1)=rand(N1,1)*pi/10;
% th(1:N1)=rand(N2,1) + 3;
% 
% phi(N1+1:end)=rand(N2,1)*pi/4+pi;
% th(N1+1:end)=rand(N2,1)*0.2 + 2;




% centered patch of density similar to steady state density
% phi=rand(N,1)*pi*2;
% th=(-rand(N,1)+1)*acosh(2);

% radially symmetric
% NN = floor(sqrt(N));
% phi = linspace(0,2*pi,NN+1);
% phi = phi(1:end-1);
% th = rand(NN,1)*2;
% [th,phi] = meshgrid(th,phi);
% th=th(:);
% phi = phi(:);

% band
% phi=rand(N,1)*2*pi;
% th=rand(N,1)*0.05 + 1;


% two symmetric patches
% N1=ceil(N/2);
% N2=N1;
% N = N1*2;
% 
% phi = zeros(N,1);
% th = zeros(N,1);
% 
% phi(1:N1)=rand(N1,1)*pi/2;
% th(1:N1)=rand(N2,1)*2;
% 
% phi(N1+1:end)=phi(1:N1)+pi;
% 
% phi=mod(phi,2*pi);
% th(N1+1:end)=th(1:N1);

% along constant phi
phi = zeros(N,1);
th = rand(N,1)*2;

%% Plotting intialization
th_up=max([th;1.4]); % upper limit of theta, determines how to plot hyperboloid
[X,Y,Z] = hyp2cart_unit(th,phi);

H = plot3(X,Y,Z,'k.'); hold on;
title('t = 0');
xy_lim = sinh(th_up)*1.2;
z_up = cosh(th_up)*1.2;
z_lo = (1-cosh(th_up)*0.1);
axis([-xy_lim, xy_lim,-xy_lim, xy_lim,z_lo,z_up]);

phi_center = mod(sum(phi),2*pi)
% th_center = acosh(mean(cosh(th)));
th_center= mean(th);
[X_center, Y_center, Z_center] = hyp2cart_unit(th_center,phi_center);
H_center = plot3(X_center, Y_center, Z_center,'r*');

%% plot hyperboloid
N_fine = 100;
[x,y,z,th_fine,phi_fine]=hyperboloid(th_up,N_fine);

% hyp = surfl(x,y,z);
% set(hyp, 'FaceAlpha', 1)
% shading interp;

hyp = surf(x,y,z);
set(hyp,'edgecolor','none');

axis square;
az = 45;
el = 20;
% el = -90;
view(az, el);
hold off;

%% plot density
figure(2);
% G = plot(th,sinh(th).*phi,'.');
[XX,YY] = pol2cart(phi,sqrt(2*(cosh(th)-1)));
[XX_center,YY_center] = pol2cart(phi_center,sqrt(2*(cosh(th_center)-1)));
G = plot(XX,YY,'.'); hold on;
G_center = plot(XX_center, YY_center,'r*');
% xlabel('\theta');
% ylabel('sin(\theta)\phi');
title('t = 0');
ax=gca;
axis(ax,'equal');

% [x_center, y_center] = pol2cart(mean(phi), mean(sqrt(2*(cosh(th)-1))));
% plot(x_center,y_center,'*');
hold off;


% figure(3);
% [XX,YY] = pol2cart(
% F =

pause;
 

% press any key to continue

%% Forward Euler
dth=zeros(N,1);
dphi=zeros(N,1);

for k=1:step_num
    
    % actual particles
    for i=1:N
        th_i=th(i);
        phi_i=phi(i);
        
        sum_th = 0;
        sum_phi = 0;
        for j=1:N
            if j==i
                continue;
            else
                
                th_j=th(j);
                phi_j=phi(j);
                
                sinh_th_i = sinh(th_i);
                sinh_th_j = sinh(th_j);
                cosh_th_i = cosh(th_i);
                cosh_th_j = cosh(th_j);
                
                z=cosh_th_i*cosh_th_j - sinh_th_i*sinh_th_j*cos(phi_i-phi_j);
                rho = acosh(z);
                
                dzdth = sinh_th_i*cosh_th_j - cosh_th_i*sinh_th_j*cos(phi_i-phi_j);
                dzdphi = sinh_th_i*sinh_th_j*sin(phi_i-phi_j);
                
                % sin potential
                % dK/dz
                %                 dKdz = -1/(4*pi)/(sinh(rho/2).*cosh(rho/2).*sqrt(z^2-1) )...
                %                     + C*coth(rho)./sqrt(z^2-1);
                
                % % variable q
                dKdz = -1/(4*pi)/(sinh(rho/2).*cosh(rho/2).*sqrt(z^2-1) )...
                    + C*coth(rho.^q).*rho.^(q-1)./sqrt(z^2-1);
                
                % % only attraction
                % dKdz = C*coth(rho.^q).*rho.^(q-1)./sqrt(z^2-1);
                
                sum_th = sum_th + dKdz.*dzdth;
                sum_phi = sum_phi + dKdz.*dzdphi./sinh_th_i;
            end
        end
        dth(i)=-sum_th/N;
        dphi(i)=-sum_phi/N;
    end
    
    
    
    % existing particles
    th = th + dt*dth;
    phi = phi + dt*dphi;
    phi=mod(phi,2*pi);
    [X,Y,Z] = hyp2cart_unit(th,phi);
    
    
    % stops if velocity is too big or small
    if norm([dth;dphi],Inf) < dt/100
        disp(['Speed very small at t = ' num2str(dt*k-1)]);
        break;
    elseif norm([dth;dphi],Inf) > 10^6
        disp(['Speed blow up at t = ' num2str(dt*k)]);
        break;
    end
    
    if mod(k,frame_skip)==0
        % plot on hyperboloid
        set(H,{'xdata','ydata','zdata'},{X,Y,Z});   
        phi_center = mean(phi);
        % th_center = acosh(mean(cosh(th)));
        th_center= mean(th)
        [X_center, Y_center, Z_center] = hyp2cart_unit(th_center,phi_center);
        set(H_center,{'xdata','ydata','zdata'},{X_center,Y_center,Z_center}); 
        title(['t = ' num2str(k*dt) ', v inf-norm = ' num2str(norm([dth;dphi],Inf))]);

        
        % particle paths
        % set(Hv,{'xdata','ydata','zdata'},{Xv,Yv,Zv});
        % line(X_path, Y_path,Z_path,'color',[1,0.5,0],'LineWidth',1,'LineStyle','-');
        
        
        % plot on density graph
        [XX,YY] = pol2cart(phi,sqrt(2*(cosh(th)-1)));
        [XX_center,YY_center] = pol2cart(phi_center,sqrt(2*(cosh(th_center)-1)));
        
        set(G,{'xdata','ydata'},{XX,YY});
        set(G_center,{'xdata','ydata'},{XX_center,YY_center});
        title(['t = ' num2str(k*dt) ', v inf-norm = ' num2str(norm([dth;dphi],Inf))]);
        
        drawnow;
    end
end

%% plot equilibirum

set(H,{'xdata','ydata','zdata'},{X,Y,Z});
title(['t = ' num2str(k*dt) ', v inf-norm = ' num2str(norm([dth;dphi],Inf))]);


