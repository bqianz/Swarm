%% 2-sphere: potential giving uniform density distribution at equilibrium

%% Radius and number of particles
clear, close all;

R=1;
N = 140;
N_vac = 0;

% pick a particle to follow its path
p_num = 7;
p = randi(N,p_num,1);


%% time-step
dt=10^-2;
tmax=60;
step_num=ceil(tmax/dt);
dt = tmax/step_num;
frame_skip = floor(step_num/100);
% frame_skip = 1;
%% initial condition and plotting it
% phi=rand(1,N)'*2*pi;
% th=rand(1,N)'*pi-(pi/2);

% phi from 0 to 2pi
phi=randn(1,N)'*2*pi;
% phi = linspace(0,2*pi,N+1); phi = (phi(1:end-1))';
% th from 0 to pi, angle from z-axis
th=randn(1,N)';
[X,Y,Z] = sph2cart(phi,pi/2-th,R); % all actual particles

% vacuum particles
phi_v = rand(1,N_vac)';
th_v = rand(1,N_vac)'*pi/2+pi/2;
[Xv,Yv,Zv] = sph2cart(phi_v,pi/2-th_v,R); % vaccum

% particle paths
th_p = th(p);
phi_p = phi(p);
path_scaling = 50;
[Xp,Yp,Zp] = sph2cart(phi_p,pi/2-th_p,R+eps*dt*path_scaling);

H = plot3(X,Y,Z,'k.'); hold on;
Hv = plot3(Xv,Yv,Zv,'r.'); % empty particles
% scatter3(Xp,Yp,Zp,'r*');
title('t = 0');
p_axis=R*1.5;
axis([-p_axis, p_axis,-p_axis, p_axis,-p_axis, p_axis]);

%% plot sphere
[x,y,z]=sphere;
Re=R-eps;
sph = surfl(x*Re, y*Re, z*Re);
set(sph, 'FaceAlpha', 1)
shading interp;
axis square;
az = 45;
el = 20;
% el = -90;
view(az, el);


%% Forward Euler
dth=zeros(N,1);
dphi=zeros(N,1);

dth_v=zeros(N_vac,1);
dphi_v=zeros(N_vac,1);
% kp=1/R^2;

% set up particle paths
X_path = zeros(p_num,frame_skip+1);
Y_path = zeros(p_num,frame_skip+1);
Z_path = zeros(p_num,frame_skip+1);
X_path(:,1) = Xp;
Y_path(:,1) = Yp;
Z_path(:,1) = Zp;

%% get first frame
set(gca,'nextplot','replacechildren','visible','off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,step_num) = 0;

%% step forward

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
                
                sin_th_i = sin(th_i);
                sin_th_j = sin(th_j);
                cos_th_i = cos(th_i);
                cos_th_j = cos(th_j);
                
                z=cos_th_i*cos_th_j + sin_th_i*sin_th_j*cos(phi_i-phi_j);
                
                % sin potential
                w = 1/(4*pi) * cot(acos(z)/2) / sqrt(1-z^2) ; % dK/dz
                % cot(angle/2) / sin(angle)
                
                
                % w = 1/(4*pi) * cot((acos(z))^2/2).*acos(z) ./ sqrt(1-z^2) ; % dK/dz
                % cot(angle/2) / sin(angle)
                
                % tan potential
                % w =  1/ (2*pi * (1-z^2) ); % simplified using double angle identity
                
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
    
    
    % empty particles
    for i=1:N_vac
        th_v_i=th_v(i);
        phi_v_i=phi_v(i);
        
        sum_th_v = 0;
        sum_phi_v = 0;
        for j=1:N
            th_j=th(j);
            phi_j=phi(j);
            
            sin_th2_i = sin(th_v_i);
            sin_th_j = sin(th_j);
            cos_th2_i = cos(th_v_i);
            cos_th_j = cos(th_j);
            
            z=cos_th2_i*cos_th_j + sin_th2_i*sin_th_j*cos(abs(phi_v_i-phi_j));
            w = 1/(4*pi) * cot(acos(z)/2) / sqrt(1-z^2) ;
            
            % sinh potential
            sum_th_v = sum_th_v + ( -sin_th2_i * cos_th_j + cos_th2_i * sin_th_j * cos(abs(phi_v_i-phi_j)) ) * w / R^2;
            sum_phi_v = sum_phi_v - sin_th2_i * sin_th_j * sin( abs (phi_v_i - phi_j)) * sign (phi_v_i - phi_j) * w / ( R^2 * sin_th2_i^2);
        end
        dth_v(i)=-sum_th/N;
        dphi_v(i)=-sum_phi/N;
    end
    
    % existing particles
    th = th + dt*dth;
    phi = phi + dt*dphi;
    [th,phi] = recaliberate(th,phi);
    [X,Y,Z] = sph2cart(phi,pi/2-th,R);
    
    % vaccum particles
    th_v = th_v + dt*dth_v;
    phi_v = phi_v + dt*phi_v;
    [th_v,phi_v] = recaliberate(th_v,phi_v);
    [Xv,Yv,Zv] = sph2cart(phi_v,pi/2-th_v,R);
    
    % characterstics paths
    th_p = th(p);
    phi_p = phi(p);
    [Xp,Yp,Zp] = sph2cart(phi_p,pi/2-th_p,R+eps*dt*path_scaling);
    
    path_ind = mod(k,frame_skip) + 1;
    X_path(:,path_ind) = Xp;
    Y_path(:,path_ind) = Yp;
    Z_path(:,path_ind) = Zp;
    
    
    
    if and(norm(dth,Inf) < 10^-2, norm(dphi,Inf) < 10^-2)
        disp(['Speed very small at t = ' num2str(dt*k-1)]);
        break;
    elseif norm([dth;dphi],Inf) > 10^6
        disp(['Speed blow up at t = ' num2str(dt*k)]);
        break;
    end
    
    if mod(k,frame_skip)==0
        set(H,{'xdata','ydata','zdata'},{X,Y,Z});
        set(Hv,{'xdata','ydata','zdata'},{Xv,Yv,Zv});
        line(X_path, Y_path,Z_path,'color',[1,0.5,0],'LineWidth',1,'LineStyle','-');
        title(['t = ' num2str(k*dt)]);
        drawnow;
        f = getframe;
        im(:,:,1,k/frame_skip) = rgb2ind(f.cdata,map,'nodither');
    end
end

imwrite(im,map,'Sphere_aggregation.gif','DelayTime',0,'LoopCount',inf)

%% plot equilibirum

% set(H,{'xdata','ydata','zdata'},{X,Y,Z});
% title(['t = ' num2str(k*dt)]);

% figure(2)
% [XX,YY] = pol2cart(phi,th);
% plot(XX,YY,'.');
% axis square, axis equal;