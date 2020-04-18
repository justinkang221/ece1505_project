function [hv] = OMP(L,z,A)
r = z;
[~,Gb] = size (A) ;
hv = (1:Gb)'.^0 - 1;
t = 0;I= [];
i=0;
while t < Gb && i < L
    [~,g] = max(abs(A'*r));
    t = t+1;I = unique([I,g]);
    x = A(:,I) \ z;
    r = z - A(:,I)*x;
    hv(I) = x;
    i=i+1;
end
hv = hv;