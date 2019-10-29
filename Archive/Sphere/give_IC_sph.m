function [phi,th] = give_IC_sph(N,IC_style)

switch IC_style
    case 'regular',
        phi=randn(1,N)'*2*pi;
        th=randn(1,N)';
    otherwise,
        phi=randn(1,N)'*2*pi;
        th=randn(1,N)';
end
        
end