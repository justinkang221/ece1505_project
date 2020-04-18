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
params.Nb = 30;
params.Gb_list = ceil(logspace(0,3,20)*params.M);
params.sigma_2 = 1;
params.Lp = 2;
params.d_lambda = 1/2;
if(params.L * params.Nrf > params.M)
    error('L * Nrf should be less than or equal to M');
end
% Signal power
SNR_dB = [10 20];
SNR = 10.^(SNR_dB/10);
P = SNR * params.sigma_2; 
%% Construct Sensing Matrix
params.W = generate_W(params);
%%
load("matlab.mat");
params.H = {};
for i=1:params.Nb
     params.H{i} = generate_chl(params.Lp, params.M, params.K, params.d_lambda);%H_set{i};
end
%% run simulation of omp (Assume Lp known)
type = "omp";
for i=1:length(params.Gb_list)
    params.Gb = params.Gb_list(i);
    mse_v_gb_omp_10dB(i) = chl_est_hybrid_grid_func(params, P(1), type);
end
%% run simulation of omp (Assume Lp known)
type = "omp";
for i=1:length(params.Gb_list)
    params.Gb = params.Gb_list(i);
    mse_v_gb_omp_20dB(i) = chl_est_hybrid_grid_func(params, P(2), type);
end
%%
atomic_mmse = [0.0299736202796685,0.0027451119793846];
mse_v_gb_atom_10dB = ones(length(params.Gb_list))*atomic_mmse(1);
mse_v_gb_atom_20dB = ones(length(params.Gb_list))*atomic_mmse(2);
%% Plots
figure
loglog(params.Gb_list, mse_v_gb_omp_20dB, "g-*", "LineWidth", 2)
hold on
loglog(params.Gb_list, mse_v_gb_atom_20dB, 'k--', "Linewidth", 2)
loglog(params.Gb_list, mse_v_gb_omp_10dB, "r-*", "LineWidth", 2)
loglog(params.Gb_list, mse_v_gb_atom_10dB, 'b--', "Linewidth", 2)
xlabel("Number of Grid Points")
ylabel("NMSE")
legend("OMP (SNR = 20dB)", "Atomic Norm (SNR = 20dB)", "OMP (10dB)", "Atomic Norm (10db)")
axis([min(params.Gb_list), max(params.Gb_list), 2e-3, 0.5])