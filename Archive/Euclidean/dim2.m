%% 2-dimensional aggregation with radial symmetric force F(r) = 1/r^(n-1) - r

clear, close all;
filename = '2D_aggregation.gif';

dim=2;
N = 100;
scale=1;
Z0 = randn(N*dim,1)*scale - scale/2;

x_center = mean(Z0(1:N));
y_center = mean(Z0(N+1:end));

range = norm([Z0(1:N) - x_center ; Z0(N+1:end) - y_center],inf);

tfinal=70;

[T_out, Z_out] = ode23(@DE_dim2, [0,tfinal], Z0);

num_steps=length(T_out);

% plot_ind=[1, round(num_steps/3), round(num_steps/3*2), num_steps];

% for i=1:4
%     subplot(2,2,i);
%     plot(Z_out(plot_ind(i),1:N), Z_out(plot_ind(i),N+1:end),'.');
%     title(['t = ' num2str(T_out(plot_ind(i)))]);
%     axis square, axis equal;
% end

% subplot(2,2,1)
% plot(Z_out(1,1:N),Z_out(end,N+1:end),'.');
% title('t = 0');
% axis equal, axis square;
%
% subplot(2,2,2)
% plot(Z_out(end,1:N),Z_out(end,N+1:end),'.');
% title(['t = ' num2str(tfinal)]);
% axis equal, axis square;

%%
framemax=length(T_out);
% F(framemax) = struct('cdata',[],'colormap',[]);
figure;
H = scatter(Z_out(1,1:N), Z_out(1,N+1:2*N),'.');
title('t =0' );
axis([x_center-range, x_center+range y_center-range y_center+range]);
axis square;
set(gca,'nextplot','replacechildren','visible','off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,20) = 0;


% drawnow;
% f = getframe;
% % F(1) = f;
% im = frame2im(f);
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif', 'Loopcount',inf);


for i=1:framemax
    set(H,{'xdata','ydata'},{Z_out(i,1:N), Z_out(i,N+1:2*N)});
    title(['t = ' num2str(T_out(i))]);
    f = getframe;
  im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    
%     % F(i)=f;
%     im = frame2im(f);
%     [imind,cm] = rgb2ind(im,256);
%     
%     % Write to the GIF File
%     imwrite(imind,cm,filename,'gif','WriteMode','append');


    
end
imwrite(im,map,'2D_aggregation.gif','DelayTime',0,'LoopCount',inf)

%% transparent sphere
% hold on;
% [x,y,z]=sphere;
% h = surfl(x*R, y*R, z*R);
% set(h, 'FaceAlpha', 0.5)
% shading interp
% hold off;

%%
% title Replay;
% movie(F,1);