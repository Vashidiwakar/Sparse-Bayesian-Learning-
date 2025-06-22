% specifying the dimensions
N = 20;
M = 40;
D0 = 7;
% w_true = randn(M,1);
w_true = zeros(M,1);
idx = randperm(M,D0);
w_true(idx) = randn(D0,1);
Phi = randn(N,M); % dictionary or design matrix

% Noise
SNR_dB = 0:5:30; % SNR -> signal to noise ratio
% SNR values from 0 to 30 dB in steps of 5
SNR_lin = 10.^(SNR_dB/10);

NMSE = zeros(size(SNR_lin));  % Store NMSE for each SNR
max_avg = 25;
for i = 1:length(SNR_lin)
    nmse_sum = 0;
    for j = 1:max_avg % Monte carlo loop
        sigma2 = 1 / SNR_lin(i);
        noise = sqrt(sigma2) * randn(N,1);
        t = signalGen(Phi, w_true, noise);
        Mn_new = SBL2(Phi, t, sigma2, M);
        w_est = Mn_new;
        nmse_temp = norm(w_est - w_true)^2 / norm(w_true)^2;
        nmse_sum = nmse_sum + nmse_temp;
    end
    NMSE(i) = nmse_sum/20;
end
% Plot NMSE vs SNR
figure;
semilogy(SNR_dB, NMSE, '-o', 'LineWidth',2);
xlabel('SNR (dB)');
ylabel('NMSE');
title('NMSE vs SNR for Sparse Bayesian Regression (EM)');
grid on;