function du = DE_dim3(t,u)
len=length(u);
N = len/3;
du=zeros(len,1);
a = 1/(4*pi);
q=3;
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
        r = sqrt((Xj- Xi)^2 + (Yj - Yi)^2 + (Zj-Zi)^2);
        sumX = sumX + (a/r^3 - r^(q-2)) * (Xi- Xj);
        sumY = sumY + (a/r^3 - r^(q-2)) * (Yi - Yj);
        sumZ = sumZ + (a/r^3 - r^(q-2)) * (Zi - Zj);
        end
    end
du(i) = sumX/N;
du(N+i) = sumY/N;
du(2*N+i) = sumZ/N;
end
end

