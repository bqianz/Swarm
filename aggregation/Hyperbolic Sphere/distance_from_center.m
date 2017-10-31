function distance_from_center(th,phi)
th_center = mean(th);
phi_center = mean(phi);

d = abs(acosh(  cosh(th)*cosh(th_center) - sinh(th)*sinh(th_center).*cos(phi-phi_center)    ));

[x,y] = pol2cart(phi,sqrt(2*(cosh(th)-1)));
x_center = mean(x);
y_center = mean(y);
% [x_center,y_center] = pol2cart(phi_center,sqrt(2*(cosh(th_center)-1)));
% 
% [x,y] = pol2cart(phi,th);
% [x_center,y_center] = pol2cart(phi_center,th_center);

figure;
plot3(x,y,d,'o'); hold on;
plot3(x_center,y_center,0,'*');


% [xq,yq] = meshgrid();
% vq = griddata(x,y,v,xq,yq);
% Plot the gridded data as a mesh and the scattered data as dots.
% 
% mesh(xq,yq,vq)
% hold on
% plot3(x,y,v,'o')
% xlim([-2.7 2.7])
% ylim([-2.7 2.7])


title('Distance from center patch');

axis equal;
end