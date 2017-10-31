function [phi,th] = give_initial(N,IC_style)

switch IC_style
    case 'regular',
        phi=rand(N,1)*2*pi;
        th=rand(N,1);
        
    case 'tall_far',
        phi=rand(N,1)*pi/10;
        th=rand(N,1) + 3;
        
    case 'tall_far_two'
        N1 = round(N/2);
        N2 = N - N1;
        
        phi = zeros(N,1);
        th = zeros(N,1);
        
        phi(1:N1)=rand(N1,1)*pi/10;
        th(1:N1)=rand(N2,1) + 3;
        
        phi(N1+1:end)=rand(N2,1)*pi/4+pi;
        th(N1+1:end)=rand(N2,1)*0.2 + 2;
        
    case 'already_uniform_density'
        phi=rand(N,1)*pi*2;
        th=(-rand(N,1)+1)*acosh(2);
        
    case 'radial_symmetric'
        NN = floor(sqrt(N));
        phi = linspace(0,2*pi,NN+1);
        phi = phi(1:end-1);
        th = rand(NN,1)*2;
        [th,phi] = meshgrid(th,phi);
        th=th(:);
        phi = phi(:);
        
    case 'band'
        phi=rand(N,1)*2*pi;
        th=rand(N,1)*0.05 + 1;
        
    case 'regular_duplicate'
        N1=ceil(N/2);
        N2=N1;
        N = N1*2;
        
        phi = zeros(N,1);
        th = zeros(N,1);
        
        phi(1:N1)=rand(N1,1)*pi*2;
        th(1:N1)=rand(N2,1)*2;
        
        phi(N1+1:end)=phi(1:N1)+pi;
        phi=mod(phi,2*pi);
        th(N1+1:end)=th(1:N1);
        
    case 'far_duplicated'
        N1 = round(N/2);
        N2 = N - N1;
        
        phi = zeros(N,1);
        th = zeros(N,1);
        
        phi(1:N1)=rand(N1,1)*pi/10;
        th(1:N1)=rand(N2,1) + 3;
        
        phi(N1+1:end)=rand(N2,1)*pi/4+pi;
        th(N1+1:end)=rand(N2,1)*0.2 + 2;
        
    otherwise, % regular
        phi=rand(N,1)*2*pi;
        th=rand(N,1);
       s 
end
end