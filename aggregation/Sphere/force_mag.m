clear;close all;
N=100;
angle_range = linspace(0,pi,N+2);
angle_range = angle_range(2:end-1);
qi=[1,0];
force_magnitude=zeros(size(angle_range));
newt_proj_magnitude=zeros(size(angle_range));
newt_magnitude=zeros(size(angle_range));

i=1;
for th=angle_range
    dt=cos(th);
    qj=[cos(th),sin(th)];
    force_magnitude(i)=norm((qj - dt*qi)/(1-dt^2)^(3/2));
    newt_proj_magnitude(i)=cos(th/2)/2/(1-cos(th));
    newt_magnitude(i)=1/2/(1-cos(th));
    i=i+1;
end

plot(angle_range, [force_magnitude;newt_proj_magnitude;newt_magnitude]);
axis([0 pi 0 10]);
legend('on the sphere','tangential projection of newtonian','newtonian');
title('Magnitude of force wrt angle of two bodies')

figure(2);
plot(angle_range, [log(log(force_magnitude./newt_proj_magnitude));...
    log(force_magnitude./newt_magnitude)]);
title('Ratio')
axis([0 pi/2 -10 1]);