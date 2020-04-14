% Performs channel estimation 
clc; clear all; close all;
% Some definitions
L = 2; K = 10; M = 64; Nrf = 16; Qb = 7;
Lp = 2; d_lambda = 1/2; Gb = 2 * M;
batch_size = 20;
P = 1;

if(L * Nrf > M)
    error('L * Nrf should be less than or equal to M');
end
% Construct Wrf randomly
c = randi(2^Qb - 1,L * Nrf, 1);
%c = zeros(L * Nrf, 1);
[X, Y] = meshgrid(0:(M-1), c);
Wrf = exp(1j * 2 * pi/2^Qb * X .* Y);
idx = randperm(M);
A = zeros(L * Nrf, M);
for ii = 1:(L * Nrf)
    A(ii, idx(ii)) = 1;
end


phase_vec = linspace(-1,1,Gb+1);
phase_vec = phase_vec(2:end);
psi = zeros(M,Gb);
for ii=1:Gb
    psi(:,ii) = exp(1i*2*pi*d_lambda*phase_vec(ii)*(0:(M-1)));
end


W = zeros(L * Nrf, L * M);
for ii = 1:L
    temp = A((ii-1)* Nrf + (1:Nrf), :).*Wrf((ii-1)* Nrf + (1:Nrf), :);
    W((ii-1)* Nrf + (1:Nrf), (ii-1)* M + (1:M)) = temp;
end

A = W * sqrt(P) * repmat(psi,[L,1]);
SNR_dB = 0:2:20;
mmse_atom = zeros(1, length(SNR_dB));
mmse_omp = zeros(1, length(SNR_dB));

for jj = 1:length(SNR_dB)
    crt_SNR = 10^(SNR_dB(jj)/10);
    sigma_2 = P/crt_SNR; % Noise power;
    mmse_atom_snr = zeros(1, batch_size);
    mmse_omp_snr = zeros(1, batch_size);
    tau = choose_best_tau(SNR_dB(jj), M, L, Lp, W);
for ii = 1:batch_size
    % Generate channel
    H = generate_chl(Lp, M, K, d_lambda);
    
    % Generate noise realization
    N = sqrt(sigma_2/2) * (randn(L * M, K) + 1j * randn(L * M, K));

    
    % Transmit pilots
    Z = W * (sqrt(P) * repmat(H,[L,1]) + N);
    
    % Channel estimation
    H_est_atom = zeros(M, K);
    H_est_omp = zeros(M, K);
    for kk = 1:K
        crt_z = Z(:, kk);
        cvx_begin sdp
            variables t u1
            variable u(M - 1) complex
            variable x(M) complex
            minimize 1/2 * sum_square_abs(sqrt(P) * W * repmat(x, [L, 1]) - crt_z) + tau/2 * (t + u1)
            subject to
                [toeplitz([u1; u]) x; x' t] >= 0
        cvx_end
        
        H_est_atom(:, kk) = x;
        H_est_omp(:, kk) = psi * OMP(Lp,crt_z, A);
    end
    mmse_atom_snr(ii) = 1/(M * K) * norm(H_est_atom - H,'fro')^2;
    mmse_omp_snr(ii) = 1/(M * K) * norm(H_est_omp - H, 'fro')^2;
end
mmse_atom(jj) = mean(mmse_atom_snr);
mmse_omp(jj) = mean(mmse_omp_snr);
end

semilogy(SNR_dB, mmse_atom, 'b-x', 'LineWidth', 3)
hold on
semilogy(SNR_dB, mmse_omp, 'g--o', 'LineWidth', 3)
legend('Atomic', 'OMP')
xlabel('SNR (dB)')
ylabel('MMSE')