function dZ = DE_dim2_test(t,Z)
len=length(Z);
N = len/2;
dZ=zeros(len,1);

for i=1:N
    Xi = Z(i);
    Yi = Z(N+i);
    sumX = 0;
    sumY = 0;
    for j=1:N
        if j==i
            continue;
        else
        Xj = Z(j);
        Yj = Z(N+j);
        r = sqrt((Xj- Xi)^2 + (Yj - Yi)^2);
        sumX = sumX + (-1) * (Xi- Xj);
        sumY = sumY + (-1) * (Yi - Yj);
        end
    end
dZ(i) = sumX/N;
dZ(N+i) = sumY/N;
end
end