function du = DE_2sphere(t,u)
len=length(u);
N = len/3;
du=zeros(len,1);
q=2;
k=1;
% k=1/(u(1)^2+u(1+N)^2+u(1+2*N)^2);
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
            kq=k*(Xi*Xj + Yi*Yj + Zi*Zj);
            attX=k*kq*(Xj-kq*Xi)/(1-kq^2);
            attY=k*kq*(Yj-kq*Yi)/(1-kq^2);
            attZ=k*kq*(Zj-kq*Zi)/(1-kq^2);
            sumX = sumX + attX/(2*pi) - 1/attX;
            sumY = sumY + attY/(2*pi) - 1/attY;
            sumZ = sumZ + attZ/(2*pi) - 1/attZ;
        end
    end
    du(i)=sumX;
    du(N+i)=sumY;
    du(2*N+i)=sumZ;
end
du=du/(N);
end

