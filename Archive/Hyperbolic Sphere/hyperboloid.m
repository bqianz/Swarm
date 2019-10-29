function [xx,yy,zz,theta,phi] = hyperboloid(varargin)
narginchk(0,2);
[cax,args,nargs] = axescheck(varargin{:});

n = 20;
theta_upper=1;
if nargs > 0
    theta_upper = args{1};
    if nargs > 1;
        n=args{2};
    end
end

% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
% theta = (-n:2:n)/n*pi;
% phi = (-n:2:n)'/n*pi/2;
theta = linspace(0,theta_upper,n+1)';
phi = linspace(0,2*pi,n+1);
sinhtheta = sinh(theta);

%[phi_mesh,th_mesh] = meshgrid(phi,theta);

x = sinhtheta*cos(phi);
y = sinhtheta*sin(phi);
z = cosh(theta)*ones(1,n+1);

if nargout == 0
    cax = newplot(cax);
    surf(x,y,z,'parent',cax)
else
    xx = x; yy = y; zz = z;
end
