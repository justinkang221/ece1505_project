function [hv] = BPD(lambda,z,A)
%CVX will also work, but this is way faster
% [~,Gb] = size (A);
% cvx_begin quiet
%     variable hv(Gb) complex
%     minimize 0.5*sum_square_abs(z - A*hv) + lambda*norm(hv,1)
% cvx_end
[~, hv] = evalc("as_bpdn(A, z, lambda)");
