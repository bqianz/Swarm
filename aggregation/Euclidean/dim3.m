%% 3-dimensional aggregation with radial symmetric force F(r) = a/r^(n-1) - r

%%
clear; close all;

R = (4*pi)^(-1/3); % steady state radius = approximately 0.4301
% a = (4*pi)^(-1)
dim=3;
N = 200;

u0 = rand(1,N*dim).*2 - 1;

[T_out, u_out] = ode23(@DE_dim3, [0,50], u0);

%% plot frames
framemax=length(T_out);
F(framemax) = struct('cdata',[],'colormap',[]);
close all;
figure;
for i=1:framemax
    scatter3(u_out(i,1:N), u_out(i,N+1:2*N),u_out(i,2*N+1:3*N),'.');
    title(['t = ' num2str(T_out(i))]);
    axis([-1 1 -1 1 -1 1]);
    F(i)=getframe;
end

%% transparent sphere
hold on;
[x,y,z]=sphere;
h = surfl(x*R, y*R, z*R); 
set(h, 'FaceAlpha', 0.5)
shading interp
hold off;

%%
% title Replay;
% movie(F,1);
