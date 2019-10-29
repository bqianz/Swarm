%% 2-sphere: potential giving uniform density distribution at equilibrium

%% Radius and number of particles
clear, close all;

R=1;
N = 140;

filename = 'obstacle_w_quiver.gif';

%% time-step
dt=10^-2;
tmax=200;
step_num=ceil(tmax/dt);
dt = tmax/step_num;
frame_skip = floor(step_num/100);
% frame_skip = 1;
%% initial condition and boundary
% phi from 0 to 2pi
phi=rand(1,N)'*2*pi;
% phi = linspace(0,2*pi,N+1); phi = (phi(1:end-1))';

% th from 0 to pi, angle from z-axis
th=rand(1,N)'*0.1;
th_bd=0.5*pi;

I = (th >= th_bd); % indices of particles that have reached the boundary

[X,Y,Z] = sph2cart(phi,pi/2-th,R);



%% Plot particles
F = figure;
H = plot3(X,Y,Z,'k.'); hold on;
title('t = 0');
p_axis=R*1.5;
axis([-p_axis, p_axis,-p_axis, p_axis,-p_axis, p_axis]);

H2 = [];

%% plot sphere
[x,y,z]=sphere;
Re=R-eps;
sph = surfl(x*Re, y*Re, z*Re);
set(sph, 'FaceAlpha', 1)
shading interp;
axis square;
az = 45;
el = 30;
% el = -90;
view(az, el);
axis off;
% set(gca,'nextplot','replacechildren','visible','off')

%% Forward Euler
dth=zeros(N,1);
dphi=zeros(N,1);
% kp=1/R^2;

for k=0:step_num
    
    % display message for velocity blow-up or reaching equilibrium
    
    
    if norm([dth;dphi],Inf) > 10^6
        disp(['Speed blow up at t = ' num2str(dt*k)]);
        break;
    elseif and( and(norm(dth,Inf) < 10^-3, norm(dphi,Inf) < 10^-3), k>0 )
        disp(['Speed very small at t = ' num2str(dt*k-1)]);
        break;
    end
    
    
    % plot
    [X,Y,Z] = sph2cart(phi,pi/2-th,R);
    
    if mod(k,frame_skip)==0 % plot frames at the frequency we want
        % plot particles
        set(H,{'xdata','ydata','zdata'},{X,Y,Z});
        title(['t = ' num2str(k*dt)]);
        
        % plot quiver
        if any(I);   % check if particles have reached the boundary
            % calculate quiver
            dX = R*sin(th).*(cos(phi+dphi)-cos(phi));
            dY = R*sin(th).*(sin(phi+dphi)-sin(phi));
            dZ = R*(cos(th+dth) - cos(th));
            
            
            % normalize quiver
            % scale = sqrt(dX.^2 + dY.^2 + dZ.^2);
            % dX = dX./scale;
            % dY = dY./scale;
            % dZ = dZ./scale;
            
            delete(H2);
            H2 = quiver3(X(I),Y(I),Z(I), dX(I), dY(I), dZ(I),'color','red');
        end
        
        drawnow;
        
        im = frame2im(getframe(F));
        [imind,cm] = rgb2ind(im,256);
        
        % Write to the GIF File
        if k == 0
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
    
    
    
    
    % calculate new dth, dphi
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
                
                % sin potential
                w = 1/(4*pi) * cot(acos(z)/2) / sqrt(1-z^2) ; % dK/dz
                % cot(angle/2) / sin(angle)
                
                
                % tan potential
                % w = - 1/ (2*pi * (1-z^2) ); % simplified using double angle identity
                
                % cos potential
                % w =  1/(4*pi) * tan(acos(z)/2) / sqrt(1-z^2);
                
                % test
                % w =  1/(4*pi) * z / sqrt(1-z^2);
                
                sum_th = sum_th + ( -sin_th_i * cos_th_j + cos_th_i * sin_th_j * cos(phi_i-phi_j) ) * w / R^2;
                
                sum_phi = sum_phi - sin_th_i * sin_th_j * sin(phi_i - phi_j) * w / ( R^2 * sin_th_i^2);
            end
        end
        dth(i)=-sum_th/N;
        dphi(i)=-sum_phi/N;
    end
    
    
    % make sure th stay in 0 to pi
    th = th + dt*dth;
    J1 = (th >=pi);
    th(J1) = 2*pi- th(J1);
    J2 = (th <=0);
    th(J2) = -th(J2);
    
    % make sure phi stay in 0 to 2pi
    phi=mod(phi+dt*dphi,2*pi);
    
    % project back to boundary, keep I value for plotting in next loop
    I = (th >= th_bd);
    th(I) = th_bd;
    
end

imwrite(im,map,'Sphere_aggregation.gif','DelayTime',0,'LoopCount',inf)

%% plot equilibirum
% 
% set(H,{'xdata','ydata','zdata'},{X,Y,Z});
% title(['t = ' num2str(k*dt)]);
