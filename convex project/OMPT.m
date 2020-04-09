function [hv] = OMPT(d,z,A)
r = z;
[~,Gb] = size (A) ;
hv = (1:Gb)'.^0 - 1;
t = 0;I= [];
while norm(r) > d && t < Gb
    [~,g] = max(abs(A'*r));
    t = t+1;I = unique([I,g]);
    x = pinv(A(:,I))*z;
    r = z - A(:,I)*x;
    hv(I) = x;
end
