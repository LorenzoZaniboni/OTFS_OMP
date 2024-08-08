clc; clear all; close all;
addpath(pwd + "/utils");
SNRdB = -10:2:20;
SNR = 10.^(SNRdB/10);
M = 32; l_max = 8; % l_max is the maximum delay 
N = 32; k_max = 4; % k_max is the maximum Doppler
Ptx = eye(M);
Prx = eye(N);
taps = 4;
iter = 20;
F_m = 1/sqrt(M)*dftmtx(M);
F_n = 1/sqrt(N)*dftmtx(N);
thrld = 3;
% User 1 position 
k_p_u1 = 1; 
l_p_u1 = 1;
% User 2 position 
k_p_u2 = 9;
l_p_u2 = 1;
% User 3 position
k_p_u3 = 1;
l_p_u3 = 5;
% User 4 position
k_p_u4 = 9;
l_p_u4 = 5;

pp_dB = 25;  % pilot power in the grid

% Vector definition
true_channel_matrix_1 = zeros(l_max, k_max);
true_channel_matrix_2 = zeros(l_max, k_max);
true_channel_matrix_3 = zeros(l_max, k_max);
true_channel_matrix_4 = zeros(l_max, k_max);
X_DD_u1 = zeros(M,N);
X_DD_u2 = zeros(M,N);
X_DD_u3 = zeros(M,N);
X_DD_u4 = zeros(M,N);

BER_MMSE_ICSI_LMMSE = zeros(1,length(SNRdB));
BER_MMSE_PCSI_LMMSE = zeros(1,length(SNRdB));

for ii = 1:iter
    % Channel parameters
    delay_taps_1 = randperm(l_max, taps)-1; % delay taps channel 1
    delay_taps_2 = randperm(l_max, taps)-1; % delay taps channel 2
    delay_taps_3 = randperm(l_max, taps)-1; % delay taps channel 3
    delay_taps_4 = randperm(l_max, taps)-1; % delay taps channel 4
    doppler_taps_1 = randperm(k_max, taps)-1; % Doppler taps channel 1
    doppler_taps_2 = randperm(k_max, taps)-1; % Doppler taps channel 2
    doppler_taps_3 = randperm(k_max, taps)-1; % Doppler taps channel 3
    doppler_taps_4 = randperm(k_max, taps)-1; % Doppler taps channel 4
    % Cyclic prefix 
    Ncp_1 = max(delay_taps_1);
    Ncp_2 = max(delay_taps_2);
    Ncp_3 = max(delay_taps_3);
    Ncp_4 = max(delay_taps_4);
    % Channel coefficients computation
    chan_coef_1 = (sqrt(1/2))*(randn(1,taps) + 1i*randn(1,taps)); % channel 1
    chan_coef_2 = (sqrt(1/2))*(randn(1,taps) + 1i*randn(1,taps)); % channel 2
    chan_coef_3 = (sqrt(1/2))*(randn(1,taps) + 1i*randn(1,taps)); % channel 3
    chan_coef_4 = (sqrt(1/2))*(randn(1,taps) + 1i*randn(1,taps)); % channel 4
    
    for t = 1:taps
        true_channel_matrix_1(delay_taps_1(t)+1, doppler_taps_1(t)+1) = chan_coef_1(t);
        true_channel_matrix_2(delay_taps_2(t)+1, doppler_taps_2(t)+1) = chan_coef_2(t);
        true_channel_matrix_3(delay_taps_3(t)+1, doppler_taps_3(t)+1) = chan_coef_3(t);
        true_channel_matrix_4(delay_taps_4(t)+1, doppler_taps_4(t)+1) = chan_coef_4(t);

    end

    % Effective channel in DD domain
    channel_eff_u1 = H_eff_creation(M,N, taps, chan_coef_1, delay_taps_1, doppler_taps_1, Prx, Ptx);
    channel_eff_u2 = H_eff_creation(M,N, taps, chan_coef_2, delay_taps_2, doppler_taps_2, Prx, Ptx);
    channel_eff_u3 = H_eff_creation(M,N, taps, chan_coef_3, delay_taps_3, doppler_taps_3, Prx, Ptx);
    channel_eff_u4 = H_eff_creation(M,N, taps, chan_coef_4, delay_taps_4, doppler_taps_4, Prx, Ptx);
    noise = sqrt(1/2)*(randn(1,M*N) + 1i*randn(1,M*N));

    for snr_iter = 1:length(SNRdB)
        disp("# of interation: " + ii +  "; SNR: " + SNRdB(snr_iter))
        X_DD_u1(k_p_u1,l_p_u1) = sqrt(10^(pp_dB/10));
        X_DD_u2(k_p_u2,l_p_u2) = sqrt(10^(pp_dB/10));
        X_DD_u3(k_p_u3,l_p_u3) = sqrt(10^(pp_dB/10));
        X_DD_u4(k_p_u4,l_p_u4) = sqrt(10^(pp_dB/10));

        % DD signal matrix form
        X_DD_u1(1:8,9:29) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],8,21)-1)+1i*(2*randi([0,1],8,21)-1));
        X_DD_u2(9:16,9:29) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],8,21)-1)+1i*(2*randi([0,1],8,21)-1));
        X_DD_u3(17:25,1:32) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],9,32)-1)+1i*(2*randi([0,1],9,32)-1));
        X_DD_u4(26:32,9:29) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],7,21)-1)+1i*(2*randi([0,1],7,21)-1));

        % TF signal matrix form
        X_TF_u1 = F_m*X_DD_u1*F_n';
        X_TF_u2 = F_m*X_DD_u2*F_n';
        X_TF_u3 = F_m*X_DD_u3*F_n';
        X_TF_u4 = F_m*X_DD_u4*F_n';
        
        % time domain matrix of transmitted signal 
        S_mat_u1 = Ptx*F_m'*X_TF_u1;
        S_mat_u2 = Ptx*F_m'*X_TF_u2;
        S_mat_u3 = Ptx*F_m'*X_TF_u3;
        S_mat_u4 = Ptx*F_m'*X_TF_u4;

        
        TxSamples_u1 = S_mat_u1(:).';
        TxSamples_u2 = S_mat_u2(:).';
        TxSamples_u3 = S_mat_u3(:).';
        TxSamples_u4 = S_mat_u4(:).';

        % CP transmittedd signal
        TxSamplesCP_u1 = [TxSamples_u1(M*N-Ncp_1+1:M*N) TxSamples_u1];
        TxSamplesCP_u2 = [TxSamples_u2(M*N-Ncp_2+1:M*N) TxSamples_u2];
        TxSamplesCP_u3 = [TxSamples_u3(M*N-Ncp_3+1:M*N) TxSamples_u3];
        TxSamplesCP_u4 = [TxSamples_u4(M*N-Ncp_4+1:M*N) TxSamples_u4];

        % Initializing Rx signal
        RxsamplesCP_u1 = 0;
        RxsamplesCP_u2 = 0;
        RxsamplesCP_u3 = 0;
        RxsamplesCP_u4 = 0;

        %% OTFS Tx symbol 
        for tap_iter = 1:taps
            RxsamplesCP_u1 = RxsamplesCP_u1 + chan_coef_1(tap_iter)*circshift(TxSamplesCP_u1.*exp(1j*2*pi/M*...
                (-Ncp_1:M*N-1)*doppler_taps_1(tap_iter)/N),[0,delay_taps_1(tap_iter)]);
            RxsamplesCP_u2 = RxsamplesCP_u2 + chan_coef_2(tap_iter)*circshift(TxSamplesCP_u2.*exp(1j*2*pi/M*...
                (-Ncp_2:M*N-1)*doppler_taps_2(tap_iter)/N),[0,delay_taps_2(tap_iter)]);
            RxsamplesCP_u3 = RxsamplesCP_u3 + chan_coef_3(tap_iter)*circshift(TxSamplesCP_u3.*exp(1j*2*pi/M*...
                (-Ncp_3:M*N-1)*doppler_taps_3(tap_iter)/N),[0,delay_taps_3(tap_iter)]);
            RxsamplesCP_u4 = RxsamplesCP_u4 + chan_coef_4(tap_iter)*circshift(TxSamplesCP_u4.*exp(1j*2*pi/M*...
                (-Ncp_4:M*N-1)*doppler_taps_4(tap_iter)/N),[0,delay_taps_4(tap_iter)]);
        end
        
        % received CP signal 
        Rxsamples = RxsamplesCP_u1(Ncp_1+1:M*N+Ncp_1) + RxsamplesCP_u2(Ncp_2+1:M*N+Ncp_2) + RxsamplesCP_u3(Ncp_3+1:M*N+Ncp_3) + RxsamplesCP_u4(Ncp_4+1:M*N+Ncp_4) + noise;
        R_mat = reshape(Rxsamples.',M,N);
        Y_TF = F_m*Prx*R_mat;
        Y_DD = F_m'*Y_TF*F_n;
        
        %% DD-domain CSI estimation for impulse pilot

        % User 1
        [h_est_u1,delay_taps_est_u1,doppler_taps_est_u1,num_est_taps_u1,estimated_channel_matrix_u1] = Pilot_impulse_CE_modified(M,N,k_p_u1,l_p_u1,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u1 = H_eff_creation(M,N,num_est_taps_u1,h_est_u1,delay_taps_est_u1,doppler_taps_est_u1,Prx,Ptx);
        
        % User 2
        [h_est_u2,delay_taps_est_u2,doppler_taps_est_u2,num_est_taps_u2,estimated_channel_matrix_u2] = Pilot_impulse_CE_modified(M,N,k_p_u2,l_p_u2,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u2 = H_eff_creation(M,N,num_est_taps_u2,h_est_u2,delay_taps_est_u2,doppler_taps_est_u2,Prx,Ptx);
        
        % User 3
        [h_est_u3,delay_taps_est_u3,doppler_taps_est_u3,num_est_taps_u3,estimated_channel_matrix_u3] = Pilot_impulse_CE_modified(M,N,k_p_u3,l_p_u3,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u3 = H_eff_creation(M,N,num_est_taps_u3,h_est_u3,delay_taps_est_u3,doppler_taps_est_u3,Prx,Ptx);

        % User 4
        [h_est_u4,delay_taps_est_u4,doppler_taps_est_u4,num_est_taps_u4,estimated_channel_matrix_u4] = Pilot_impulse_CE_modified(M,N,k_p_u4,l_p_u4,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u4 = H_eff_creation(M,N,num_est_taps_u4,h_est_u4,delay_taps_est_u4,doppler_taps_est_u4,Prx,Ptx);

        %% Embedded Data Detection
        x_DD_u1 = X_DD_u1(:); % reshape of the matrix into column vector
        x_DD_u2 = X_DD_u2(:); % reshape of the matrix into column vector
        x_DD_u3 = X_DD_u3(:); % reshape of the matrix into column vector
        x_DD_u4 = X_DD_u4(:); % reshape of the matrix into column vector


        x_DD_comp = [x_DD_u1;x_DD_u2;x_DD_u3;x_DD_u4]; % entire column vector

        set_S_comp = find(x_DD_comp);

        ind_u1 = find(set_S_comp == 1);
        set_S_comp(ind_u1) = [];
        ind_u2 = find(set_S_comp == 1033);
        set_S_comp(ind_u2) = [];
        ind_u3 = find(set_S_comp == 2177);
        set_S_comp(ind_u3) = [];
        ind_u4 = find(set_S_comp == 3209);
        set_S_comp(ind_u4) = [];
        % generate a vector that is composed only by data
        x_DD_S_comp = x_DD_comp(set_S_comp); 

        y_DD_vector = Y_DD(:);        
        
        %% Estimated CSI results
        y_DD_S_comp = y_DD_vector - (sqrt(10^(pp_dB/10))*H_eff_est_u1(:,1)) - (sqrt(10^(pp_dB/10))*H_eff_est_u2(:,9)) - (sqrt(10^(pp_dB/10))*H_eff_est_u3(:,129)) - (sqrt(10^(pp_dB/10))*H_eff_est_u4(:,137));
        H_eff_est_comp = [H_eff_est_u1 H_eff_est_u2 H_eff_est_u3 H_eff_est_u4];
        H_eff_est_S_comp = H_eff_est_comp(:,set_S_comp);
        x_DD_S_hat_LMMSE_comp = inv(H_eff_est_S_comp'*H_eff_est_S_comp + eye(length(set_S_comp))/SNR(snr_iter))*H_eff_est_S_comp'*y_DD_S_comp;
        DecodedBits_LMMSE_real = (real(x_DD_S_hat_LMMSE_comp)>=0);
        DecodedBits_LMMSE_imag = (imag(x_DD_S_hat_LMMSE_comp)>=0);
        DecodedSyms_LMMSE_real = sqrt((1/2)*SNR(snr_iter))*(2*DecodedBits_LMMSE_real-1);
        DecodedSyms_LMMSE_imag = sqrt((1/2)*SNR(snr_iter))*(2*DecodedBits_LMMSE_imag-1);
        BER_MMSE_ICSI_LMMSE(snr_iter) = BER_MMSE_ICSI_LMMSE(snr_iter) + sum(DecodedSyms_LMMSE_real ~= real(x_DD_S_comp))+ sum(DecodedSyms_LMMSE_imag ~= imag(x_DD_S_comp));



        %% Perfect CSI results
        y_DD_S_comp = y_DD_vector - (sqrt(10^(pp_dB/10))*channel_eff_u1(:,1)) - (sqrt(10^(pp_dB/10))*channel_eff_u2(:,9)) - (sqrt(10^(pp_dB/10))*channel_eff_u3(:,129)) - (sqrt(10^(pp_dB/10))*channel_eff_u4(:,137));   
        H_eff_comp = [channel_eff_u1 channel_eff_u2 channel_eff_u3 channel_eff_u4];       
        H_eff_S_comp = H_eff_comp(:,set_S_comp);
        x_DD_S_hat_LMMSE_comp = inv(H_eff_S_comp'*H_eff_S_comp + eye(length(set_S_comp))/SNR(snr_iter))*H_eff_S_comp'*y_DD_S_comp;
        DecodedBits_LMMSE_real = (real(x_DD_S_hat_LMMSE_comp) >= 0);
        DecodedBits_LMMSE_imag = (imag(x_DD_S_hat_LMMSE_comp) >= 0);
        DecodedSyms_LMMSE_real = sqrt((1/2)*SNR(snr_iter))*(2*DecodedBits_LMMSE_real - 1);
        DecodedSyms_LMMSE_imag = sqrt((1/2)*SNR(snr_iter))*(2*DecodedBits_LMMSE_imag - 1);
        BER_MMSE_PCSI_LMMSE(snr_iter) = BER_MMSE_PCSI_LMMSE(snr_iter) + sum(DecodedSyms_LMMSE_real ~= real(x_DD_S_comp))+ sum(DecodedSyms_LMMSE_imag ~= imag(x_DD_S_comp));

    end

end

BER_MMSE_ICSI_LMMSE = BER_MMSE_ICSI_LMMSE/(length(set_S_comp)*iter*2);
BER_MMSE_PCSI_LMMSE = BER_MMSE_PCSI_LMMSE/(length(set_S_comp)*iter*2);

figure;
semilogy(SNRdB,BER_MMSE_ICSI_LMMSE,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',4.0);
hold on;
semilogy(SNRdB,BER_MMSE_PCSI_LMMSE,'r-o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',4.0);
hold on;
grid on;
%axis tight;
%title('LMMSE');
xlabel('SNR [dB]');
ylabel('BER');
legend('Imperfect CSI','Perfect CSI')
