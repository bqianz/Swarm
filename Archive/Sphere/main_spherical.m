%% 2-sphere: potential giving uniform density distribution at equilibrium

%% Radius and number of particles
clear, close all;

R=1;
N = 1000;

% pick a particle to follow its path
p_num = 7;
p = randi(N,p_num,1);
tic;

%% time-step
dt=10^-2;
tmax=300;
step_num=ceil(tmax/dt);
dt = tmax/step_num;
frame_skip = floor(step_num/10000);
% frame_skip = 1;
%% initial condition and plotting it
[phi,th] = give_IC_sph(N,'regular');
[X,Y,Z] = sph2cart(phi,pi/2-th,R);

% density = colour_density_sph(th,phi);
% H = scatter3(X,Y,Z,5,density*50,'filled'); hold on;
% title('t = 0');
% p_axis=R*1.5;
% axis([-p_axis, p_axis,-p_axis, p_axis,-p_axis, p_axis]);



%% plot sphere
[x,y,z,thh,phh]=sphere_grid;
Re=R-eps;
sph = surfl(x*Re, y*Re, z*Re);
set(sph, 'FaceAlpha', 0.5)
shading interp;
axis square;
az = 45;
el = 20;
% el = -90;
view(az, el);

pause;

%% Forward Euler
dth=zeros(N,1);
dphi=zeros(N,1);


%% get first frame
% set(gca,'nextplot','replacechildren','visible','off')
% f = getframe;
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% % im(1,1,1,step_num) = 0;

%% step forward
for k=1:step_num
    
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
                
                sin_th_i = sin(th_i);
                sin_th_j = sin(th_j);
                cos_th_i = cos(th_i);
                cos_th_j = cos(th_j);
                
                z=cos_th_i*cos_th_j + sin_th_i*sin_th_j*cos(phi_i-phi_j);
                
                w = 1/(4*pi) * cot(acos(z)/2) / sqrt(1-z^2) ; % dK/dz

                
                sum_th = sum_th + ( -sin_th_i * cos_th_j + cos_th_i * sin_th_j * cos(phi_i-phi_j) ) * w / R^2;
                sum_phi = sum_phi - sin_th_i * sin_th_j * sin(phi_i - phi_j) * w / ( R^2 * sin_th_i^2);
            end
        end
        dth(i)=-sum_th/N;
        dphi(i)=-sum_phi/N;
    end
    
    % existing particles
    th = th + dt*dth;
    phi = phi + dt*dphi;
    [th,phi] = recaliberate(th,phi);
    [X,Y,Z] = sph2cart(phi,pi/2-th,R);
    
    
    
    
    if norm([dth;dphi],inf) < dt/10
        disp(['Speed very small at t = ' num2str(dt*k-1)]);
        break;
    elseif norm([dth;dphi],Inf) > 10^6
        disp(['Speed blow up at t = ' num2str(dt*k)]);
        break;
    end
    
    if mod(k,frame_skip)==0
        k*dt
        toc
%         set(H,{'xdata','ydata','zdata'},{X,Y,Z});
%         line(X_path, Y_path,Z_path,'color',[1,0.5,0],'LineWidth',1,'LineStyle','-');
%         title(['t = ' num2str(k*dt)]);
%         drawnow;
%         f = getframe;
%         im(:,:,1,k/frame_skip) = rgb2ind(f.cdata,map,'nodither');
    end
end

% imwrite(im,map,'Sphere_aggregation.gif','DelayTime',0,'LoopCount',inf)

%% plot equilibirum



% set(H,{'xdata','ydata','zdata'},{X,Y,Z});
% title(['t = ' num2str(k*dt)]);
% 
% % figure(2)
% % [XX,YY] = pol2cart(phi,th);
% % plot(XX,YY,'.');
% % axis square, axis equal;
% 
% %% colour density
% figure(1);
% density = colour_density_sph(th,phi);
% H = scatter3(X,Y,Z,5,density*50,'filled');
% 


%%
[mean_short, nearest_dist,farthest_dist] = mean_shortest_distance(th,phi);
density = colour_density_sph(th,phi,mean_short);

hold on;
H = scatter3(X,Y,Z,5,density*50,'filled');
title(['t = ' num2str(k*dt)]);

%%
figure(2);
plot(1:N,density,'.');
title('Density at equilibrium');
ylabel('density');

figure(3)
plot(1:N, nearest_dist,'.');
title('Distance to nearest particle');
ylabel('distance');




