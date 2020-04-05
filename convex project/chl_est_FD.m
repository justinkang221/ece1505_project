% Performs channel estimation 
clc; clear all; close all;

% Some definitions
L = 10; K = 10; M = 64; %L >= K
Lp = 2; d_lambda = 1/2;
batch_size = 20;
SNR_dB = 20;
SNR = 10^(SNR_dB/10);
sigma_2 = 1; % Noise power
P = SNR * sigma_2; % Signal power


dbg = 0;

% For simplicity, assume the rows of X are taken from the rows of eye(L)
IL = eye(L);
phi = IL(1:K, :);
X = sqrt(P/K) * phi;
mms_err = zeros(1, batch_size);
if(dbg)
    tau = 0;
else
    tau = sqrt(sigma_2) * ( 1 + 1/log(M)) * sqrt(M * log(M) + M * log(4 * pi * log(M))); 
end

for ii = 1:batch_size
    % Generate channel
    H = generate_chl(Lp, M, K, d_lambda);
    
    % Generate noise realization
    if(dbg)
        N = 0;
    else
        N = sqrt(sigma_2/2) * (randn(M, L) + 1j * randn(M, L));
    end
    
    % Transmit pilots
    Y = H * X + N;
    
    % Channel estimation
    Z = Y * phi';
    
    H_est = zeros(M, K);
    for k = 1:K
        crt_z = Z(:, k);
        cvx_begin sdp quiet
            variables t u1
            variable u(M - 1) complex
            variable x(M) complex
            minimize sum_square_abs(sqrt(P/K) * x - crt_z) + tau * (t + u1)
            subject to
                [toeplitz([u1; u]) x; x' t] >= 0
        cvx_end
        H_est(:, k) = x;
    end
    mms_err(ii) = 1/(M * K) * norm(H_est - H,'fro')^2;
    if(dbg)
        keyboard;
    end
end