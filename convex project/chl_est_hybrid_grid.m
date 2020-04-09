%% Performs channel estimation using A on-the-grid approach.
% This script is the master script for the OMP method of Channel Estimation
% 
clc; clear; close all;
%% Define the parameter struct
% The value of this is that we can save it and reproduce results.
params.L = 2;
params.K = 10;
params.M = 64;
params.Nrf = 16;
params.Qb = 7;
params.dbg = 0;
params.Nb = 100;
params.Gb = 4*params.M;
params.sigma_2 = 1;
params.Lp = 2;
params.d_lambda = 1/2;
if(params.L * params.Nrf > params.M)
    error('L * Nrf should be less than or equal to M');
end
SNR_dB = 20:-2:0;
SNR = 10.^(SNR_dB/10);
P = SNR * params.sigma_2; % Signal power
%% Construct Sensing Matrix
params.W = generate_W(params);
%% Run simulation of bpd FIX-BPD should be better than OMP
type = "bpd";
lambda_space = logspace(-1,2,15);
best_lambdas = [10,10,10,10,9,9,9,9,13,14,15];
params.lambda_list = lambda_space(best_lambdas);
for i_snr = 1:length(SNR)
    params.curr_lambda = params.lambda_list(i_snr);
    mse_v_snr_bpd(i_snr) = chl_est_hybrid_grid_func(params,P(i_snr),type);
    
end
%% run simulation of ompt
type = "ompt";
dspace = linspace(1,10,30); best_d = [18,18,17,18,17,18,18,19,20,25,28];
param.thresh_list = dspace(best_d);
for i_snr = 1:length(SNR)
    params.curr_thresh = param.thresh_list(i_snr);
    mse_v_snr_ompt(i_snr) = chl_est_hybrid_grid_func(params, P(i_snr), type);
end
%% run simulation of omp (Assume Lp known)
type = "omp";
for i_snr=1:length(SNR)
    mse_v_snr_omp(i_snr) = chl_est_hybrid_grid_func(params, P(i_snr), type);
end
%% Plots
figure
semilogy(SNR_dB, mse_v_snr_ompt, "o-");
hold on
semilogy(SNR_dB, mse_v_snr_omp, "*-")
semilogy(SNR_dB, mse_v_snr_bpd, "^-")
xlabel("SNR (dB)")
ylabel("MSE")
legend("Threshold OMP", "K-Sparsity OMP", "BPD")