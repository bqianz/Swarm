function res=hyp_patch_area(th1,th2,phi1,phi2)
% res = abs(phi2-phi1)./12.*abs( (1+4* (sinh(th1)).^2 ).^(3/2)...
%     -  (1+4* (sinh(th2)).^2 ).^(3/2) ) ;
res = abs(phi2-phi1).* abs(cosh(th2) - cosh(th1));
end