function value = C1(u0,phi)

[M,N,color] = size(phi);

numerator   = 0;
denominator = 0;

for(i=1:M)
    for(j=1:N)
        if(phi(i,j) < 0)
            numerator   = numerator + u0(i,j);
            denominator = denominator + 1;
        end
    end
end
denominator = denominator;
value = numerator / denominator;