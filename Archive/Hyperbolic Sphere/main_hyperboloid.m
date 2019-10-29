%% 2-sphere: potential giving uniform density distribution at equilibrium
clear, close all;
% global N th_fine phi_fine hyp th_coarse phi_coarse Area n N_fine;

%% Radius and number of particles

tic;
N = 200;
C = 1/(4*pi);
q = 1;

%% time-step
dt=10^-1; % initial dt
tmax=10^8;
step_num=ceil(tmax/dt);
% dt = tmax/step_num;
% frame_skip = floor(step_num/200);
frame_skip = 100;
%% initial condition

[phi,th] = give_initial(N,'tall_far');


%% Plotting intial state
th_up=max([th;1.4]); % upper limit of theta, determines how to plot hyperboloid
[X,Y,Z] = hyp2cart_unit(th,phi);

H = plot3(X,Y,Z,'k.'); hold on;
title('t = 0');
xy_lim = sinh(th_up)*1.2;
z_up = cosh(th_up)*1.2;
z_lo = (1-cosh(th_up)*0.1);
axis([-xy_lim, xy_lim,-xy_lim, xy_lim,z_lo,z_up]);

%% plot hyperboloid

hold off;

pause;



%% Forward Euler
dth=zeros(N,1);
dphi=zeros(N,1);

t=0;

v_norm = 1;

k = 1; %step number


while t < tmax;
    if v_norm < 10^-3
        dt = 0.2;
    end
    
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
                
                dKdz = -1/(4*pi)/(sinh(rho/2).*cosh(rho/2).*sqrt(z^2-1) )...
                    + C*coth(rho.^q).*rho.^(q-1)./sqrt(z^2-1);
                
                sum_th = sum_th + dKdz.*dzdth;
                sum_phi = sum_phi + dKdz.*dzdphi./sinh_th_i;
            end
        end
        dth(i)=-sum_th/N;
        dphi(i)=-sum_phi/N;
    end
    
    k=k+1;
    t = t + dt;
    th = th + dt*dth;
    phi = phi + dt*dphi;
    phi=mod(phi,2*pi);
    [X,Y,Z] = hyp2cart_unit(th,phi);
    
    if mod(k,10)==0
        v_norm = norm([dth;dphi]);
        disp(['v_inf = ' num2str(v_norm)]);
        disp(['t = ' num2str(t)]);
        disp([num2str(k) '-th step']);
        disp(['dt = ' num2str(dt)]);
        toc
        % stops if velocity is too big or small
        if v_norm < 10^-4;
            disp(['Speed very small at t = ' num2str(dt*k-1)]);
            break;
        elseif v_norm > 10^6
            disp(['Speed blow up at t = ' num2str(dt*k)]);
            break;
        end
    end
    
end

%% plot equilibirum
[mean_short, nearest_dist,~] = mean_shortest_distance_hyp(th,phi);
density = colour_density_hyp(th,phi,mean_short);


figure(1);
[x,y,z,th_fine,phi_fine]=hyperboloid(th_up);

hyp = surfl(x,y,z);
% set(hyp, 'FaceAlpha', 1)
shading interp;
set(hyp,'edgecolor','none');

axis square;

az = 45;
el = 20;
% el = -90;
view(az, el);
hold on;
H = scatter3(X,Y,Z,5,density*1000,'filled');
title(['t = ' num2str(t) ', v inf-norm = ' num2str(norm([dth;dphi],Inf))]);


figure;
plot(1:N,density,'.');
title('Density at equilibrium');
ylabel('density');




