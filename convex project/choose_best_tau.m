function tau_best = choose_best_tau(SNR_dB, M, L, Lp, W)
filename = sprintf('tau/best_tau_%d_%d_%d_%d.mat', SNR_dB, M, L, Lp);
if(isfile(filename))
    load(filename, 'tau_best')
    return;
else
K = 2;
SNR = 10^(SNR_dB/10);
P = 1;
sigma_2 = P/SNR;
tau_start = 0.1;
tau_end = 20;
tau_step = 1;
d_lambda = 0.5;

iter = 0;
best_obj = inf;
while(iter < 3)
    
for tau = tau_start:tau_step:tau_end
    H = generate_chl(Lp, M, 2, d_lambda);
    N = sqrt(sigma_2/2) * (randn(L * M, K) + 1j * randn(L * M, K));
    Z = W * (sqrt(P) * repmat(H,[L,1]) + N);
    H_est = zeros(M, 2);
    for ii = 1:K
       crt_z = Z(:, ii);
        cvx_begin sdp
            variables t u1
            variable u(M - 1) complex
            variable x(M) complex
            minimize 1/2 * sum_square_abs(sqrt(P) * W * repmat(x, [L, 1]) - crt_z) + tau/2 * (t + u1)
            %minimize (t + u1)
            subject to
                [toeplitz([u1; u]) x; x' t] >= 0
        cvx_end
        H_est(:, ii) = x;
    end
    err_crt = norm(H - H_est,'fro')^2;
    if(err_crt < best_obj)
        best_obj = err_crt;
        tau_best = tau;
    end 
end

tau_step = tau_step/10;
tau_start = max(tau_best - 20 * tau_step, 0.01);
tau_end = tau_best + 20 * tau_step;
iter = iter + 1;
end
    
    save(filename, 'tau_best')
end

