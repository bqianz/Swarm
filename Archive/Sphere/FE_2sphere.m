function [T_out, u_out] = FE_2sphere(tmax,dt,u,th_lim,R)
len=length(u);
N = len/3;
du=zeros(1,len);
% kp=1;
% k=1/(u(1)^2+u(1+N)^2+u(1+2*N)^2);
tStepN=floor(tmax/dt)+1;
T_out=zeros(tStepN,1);
u_out=zeros(tStepN,len);
t=0;
for k=1:tStepN-1
    T_out(k)=t;
    u_out(k,:)=u;
    for i=1:N
        Xi = u(i);
        Yi = u(N+i);
        Zi = u(2*N+i);
        sumX = 0;
        sumY = 0;
        sumZ = 0;
        for j=1:N
            if j==i
                continue;
            else
                Xj = u(j);
                Yj = u(N+j);
                Zj = u(2*N+j);
                
                [~,th,~]=cart2sph(Xj,Yj,Zj);
                
                if th >= th_lim
                    sumX = 0;
                    sumY = 0;
                    sumZ = 0;
                else
                    qq=(Xi*Xj + Yi*Yj + Zi*Zj);
                    % sin potential
%                     sumX = sumX + (qq*Xi - Xj)/(1-qq);
%                     sumY = sumY + (qq*Yi - Yj)/(1-qq);
%                     sumZ = sumZ + (qq*Zi - Zj)/(1-qq);
                    
                    
                    % cot potential
%                     sumX = sumX + qq*(qq*Xi - Xj)/(1-qq^2);
%                     sumY = sumY + qq*(qq*Yi - Yj)/(1-qq^2);
%                     sumZ = sumZ + qq*(qq*Zi - Zj)/(1-qq^2);

                    sumX = sumX + (qq*Xi - Xj)/(1-qq^2);
                    sumY = sumY + (qq*Yi - Yj)/(1-qq^2);
                    sumZ = sumZ + (qq*Zi - Zj)/(1-qq^2);
                end
            end
        end
        du(i)=sumX;
        du(N+i)=sumY;
        du(2*N+i)=sumZ;
    end
    du=du/(2*pi*N);
    u=u+dt*du;
    % [u(1:N),u(N+1:2*N),u(2*N+1:end), ~] = cpSphere(u(1:N),u(N+1:2*N),u(2*N+1:end),1 , [0,0,0]);
    
    [GA,TH,~]=cart2sph(u(1:N),u(N+1:2*N),u(2*N+1:end));
    
    ind = TH > th_lim;
    TH(ind) = th_lim;
    
    [u(1:N),u(N+1:2*N),u(2*N+1:end)]=sph2cart(GA,TH,R);
    
    t=t+dt;
    
%     if norm(du,Inf) < 10^-5
%         T_out=T_out(k+1,1);
%         u_out=u_out(k+1,len);
%         break
%     end
end

T_out(end)=t;
u_out(end,:)=u;
end

