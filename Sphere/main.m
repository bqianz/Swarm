%% Attraction Repulsion on 2-sphere

%%
clear, close all;

th_lim = pi/2 + eps; % no hole
% th_lim=pi/4; % with hole

dt=10^-1;
R=1;
N = 100;

% [X,Y,Z] = sph2cart(rand(1,N)*pi,rand(1,N)*pi/2-(pi/2),R);
[X,Y,Z] = sph2cart(rand(1,N)*2*pi,pi/2-rand(1,N)*0.2,R);
u=[X,Y,Z];

%%

% [T_out, u_out] = ode23(@DE_2sphere, [0,0.02], u);

[T_out, u_out] = FE_2sphere(25,dt, u, th_lim, R);

%% plot
figure; hold on;
H=plot3(X,Y,Z,'.');
title('t = 0');
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);

[x,y,z]=sphere;
Re=R-eps;
h = surfl(x*Re, y*Re, z*Re);
set(h, 'FaceAlpha', 1)
shading interp;
axis square;
az = 45;
el = 20;
view(az, el);

framemax=length(T_out);
F(framemax) = struct('cdata',[],'colormap',[]);

for i=1:framemax
    set(H,{'xdata','ydata','zdata'},{u_out(i,1:N), u_out(i,N+1:2*N),u_out(i,2*N+1:3*N)});
    title(['t = ' num2str(T_out(i))]);
    drawnow;
    F(i)=getframe;
end


%%
% title('Cotangent Potential');
% movie(F,1);
