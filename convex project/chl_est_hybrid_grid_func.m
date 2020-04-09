function avg_mse = chl_est_hybrid_grid_func(params, P, type)
Gb = params.Gb;
Nb = params.Nb;
M = params.M;
K = params.K;
dbg = params.dbg;
sigma_2 = params.sigma_2;
L = params.L;
W = params.W;
Lp = params.Lp;
d_lambda = params.d_lambda;
if type=="ompt"
    d = params.curr_thresh;
elseif type=="bpd"
    lambda = params.curr_lambda;
end
mms_err = zeros(1, Nb);
phase_vec = linspace(-1,1,Gb+1);
phase_vec = phase_vec(2:end);
psi = zeros(M,Gb);
for i=1:Gb
    psi(:,i) = exp(1i*2*pi*d_lambda*phase_vec(i)*(0:(M-1)));
end
for ii = 1:Nb
    % Generate channel
    H = generate_chl(Lp, M, K, d_lambda);
    
    % Generate noise realization
    if(dbg)
        N = 0;
    else
        N = sqrt(sigma_2/2) * (randn(L * M, K) + 1j * randn(L * M, K));
    end
    
    % Transmit pilots
    Z = W * (sqrt(P/(K * L)) * repmat(H,[L,1]) + N);
    
    % Channel estimation
    H_est = zeros(M, K);
    for k = 1:K
        crt_z = Z(:, k);
        A = W * sqrt(P/(K * L)) * repmat(psi,[L,1]);
        if type == "ompt"
            hv = OMPT(d, crt_z, A);
        elseif type == "omp"
            hv = OMP(Lp, crt_z, A);
        elseif type == "bpd"
            hv = BPD(lambda, crt_z, A);
        else
            keyboard;
        end
        H_est(:, k) = psi*hv;
    end
    mms_err(ii) = 1/(M * K) * norm(H_est - H,'fro')^2;
end
avg_mse = sum(mms_err)/Nb;
