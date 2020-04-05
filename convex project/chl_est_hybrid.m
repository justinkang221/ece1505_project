% Performs channel estimation 
clc; clear all; close all;

% Some definitions
L = 20; K = 10; M = 64; Nrf = 2; Qb = 7;
Lp = 2; d_lambda = 1/2;
batch_size = 1e3;
SNR_dB = 20;
SNR = 10^(SNR_dB/10);
sigma_2 = 1; % Noise power
P = SNR * sigma_2; % Signal power

if(L * Nrf > M)
    error('L * Nrf should be less than or equal to M');
end
% Construct Wrf randomly
c = randi(2^Qb - 1,L * Nrf, 1);
[X, Y] = meshgrid(0:(M-1), c);
Wrf = exp(1j * 2 * pi/2^Qb * X .* Y);
idx = randperm(M);
A = zeros(L * Nrf, M);
for ii = 1:(L * Nrf)
    A(ii, idx(ii)) = 1;
end

W = zeros(L * Nrf, L * M);
for ii = 1:L
    temp = A((ii-1)* Nrf + (1:Nrf), :).*Wrf((ii-1)* Nrf + (1:Nrf), :);
    W((ii-1)* Nrf + (1:Nrf), (ii-1)* M + (1:M)) = temp;
end

dbg = 0;
fact = 1;
if(dbg)
    tau = 0;
else
    tau = sqrt(sigma_2) * ( 1 + 1/log(Nrf * L)) * sqrt((Nrf * L) * log(Nrf * L) + (Nrf * L) * log(4 * pi * log(Nrf * L))); 
end

%if(dbg)
%    eta = 1e-4;
%else
%    eta = L * Nrf * sigma_2 + fact * sqrt(L * Nrf) * sigma_2;
%end
% For simplicity, assume the rows of X are taken from the rows of eye(L)
mms_err = zeros(1, batch_size);

for ii = 1:batch_size
    % Generate channel
    H = generate_chl(Lp, M, K, d_lambda);
    
    % Generate noise realization
    if(dbg)
        N = 0;
    else
        N = sqrt(sigma_2/2) * (randn(L * M, K) + 1j * randn(L * M, K));
    end
    
    % Transmit pilots
    Z = W * (repmat(H,[L,1]) + N);
    
    % Channel estimation
    H_est = zeros(M, K);
    for k = 1:K
        crt_z = Z(:, k);
        cvx_begin sdp
            variables t u1
            variable u(M - 1) complex
            variable x(M) complex
            minimize sum_square_abs(W * repmat(x,[L,1]) - crt_z) + tau * (t + u1)
            subject to
                [toeplitz([u1; u]) x; x' t] >= 0
                %sum_square_abs(W * repmat(x,[L,1]) - crt_z) <= eta
        cvx_end
        H_est(:, k) = x;
    end
    mms_err(ii) = 1/(M * K) * norm(H_est - H,'fro')^2;
    if(dbg)
        keyboard;
    end
end