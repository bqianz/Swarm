function [TH,PHI] = recaliberate(th,phi)
    J1 = (th >=pi);
    th(J1) = 2*pi- th(J1);
    J2 = (th <=0);
    th(J2) = -th(J2);
    TH = th;
    PHI=mod(phi,2*pi);
end