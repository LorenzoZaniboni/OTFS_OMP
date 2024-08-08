

% FOUR USER CASE

clc; clear all; close all;
SNRdB = -10:1:20;
% SNRdB = 10;
SNR = 10.^(SNRdB/10);
M = 32; l_max = 8;
N = 32; k_max = 4;
Ptx = eye(M);
Prx = eye(N);
taps = 4;
iter = 20;
F_m = 1/sqrt(M)*dftmtx(M);
F_n = 1/sqrt(N)*dftmtx(N);
thrld = 3;
k_p_u1 = 1;
l_p_u1 = 1;
k_p_u2 = 9;
l_p_u2 = 1;
k_p_u3 = 1;
l_p_u3 = 5;
k_p_u4 = 9;
l_p_u4 = 5;
pp_dB = 25;  %pilot power in the grid...
BER_MMSE_ICSI_LMMSE = zeros(1,length(SNRdB));
BER_MMSE_PCSI_LMMSE = zeros(1,length(SNRdB));

% mse_impulse_pilot_u1 = zeros(1,length(SNR));
% mse_impulse_pilot_u2 = zeros(1,length(SNR));
% mse_impulse_pilot_u3 = zeros(1,length(SNR));
% mse_impulse_pilot_u4 = zeros(1,length(SNR));



for n = 1:iter
    n
    delay_taps_1 = randperm(l_max,taps)-1;
    doppler_taps_1 = randperm(k_max,taps)-1;
    delay_taps_2 = randperm(l_max,taps)-1;
    doppler_taps_2 = randperm(k_max,taps)-1;
    delay_taps_3 = randperm(l_max,taps)-1;
    doppler_taps_3 = randperm(k_max,taps)-1;
    delay_taps_4 = randperm(l_max,taps)-1;
    doppler_taps_4 = randperm(k_max,taps)-1;
    Ncp_1 = max(delay_taps_1);
    Ncp_2 = max(delay_taps_2);
    Ncp_3 = max(delay_taps_3);
    Ncp_4 = max(delay_taps_4);
    chan_coef_1 = sqrt(1/2)*(randn(1,taps)+ 1i*randn(1,taps));
    chan_coef_2 = sqrt(1/2)*(randn(1,taps)+ 1i*randn(1,taps));
    chan_coef_3 = sqrt(1/2)*(randn(1,taps)+ 1i*randn(1,taps));
    chan_coef_4 = sqrt(1/2)*(randn(1,taps)+ 1i*randn(1,taps));
    true_channel_matrix_1 = zeros(l_max,k_max);
    true_channel_matrix_2 = zeros(l_max,k_max);
    true_channel_matrix_3 = zeros(l_max,k_max);
    true_channel_matrix_4 = zeros(l_max,k_max);
    for t = 1:taps
        true_channel_matrix_1(delay_taps_1(t)+1,doppler_taps_1(t)+1) = chan_coef_1(t);  %since zero tap is considered, we need (0,0) possibility and is achieved by incrementing
        true_channel_matrix_2(delay_taps_2(t)+1,doppler_taps_2(t)+1) = chan_coef_2(t);
        true_channel_matrix_3(delay_taps_3(t)+1,doppler_taps_3(t)+1) = chan_coef_3(t);
        true_channel_matrix_4(delay_taps_4(t)+1,doppler_taps_4(t)+1) = chan_coef_4(t);
    end

    H_eff_u1 = H_eff_creation(M,N,taps,chan_coef_1,delay_taps_1,doppler_taps_1,Prx,Ptx);  %delay doppler domain actual channel 
    H_eff_u2 = H_eff_creation(M,N,taps,chan_coef_2,delay_taps_2,doppler_taps_2,Prx,Ptx);
    H_eff_u3 = H_eff_creation(M,N,taps,chan_coef_3,delay_taps_3,doppler_taps_3,Prx,Ptx);
    H_eff_u4 = H_eff_creation(M,N,taps,chan_coef_4,delay_taps_4,doppler_taps_4,Prx,Ptx);
    noise = sqrt(1/2)*(randn(1,M*N) + 1i*randn(1,M*N));    
    for snr_iter = 1:length(SNR)
        X_DD_u1 = zeros(M,N);
        X_DD_u2 = zeros(M,N);
        X_DD_u3 = zeros(M,N);
        X_DD_u4 = zeros(M,N);

        X_DD_u1(k_p_u1,l_p_u1) = sqrt(10^(pp_dB/10));
        X_DD_u2(k_p_u2,l_p_u2) = sqrt(10^(pp_dB/10));
        X_DD_u3(k_p_u3,l_p_u3) = sqrt(10^(pp_dB/10));
        X_DD_u4(k_p_u4,l_p_u4) = sqrt(10^(pp_dB/10));
       
        X_DD_u1(1:8,9:29) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],8,21)-1)+1i*(2*randi([0,1],8,21)-1));
        X_DD_u2(9:16,9:29) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],8,21)-1)+1i*(2*randi([0,1],8,21)-1));
        X_DD_u3(17:25,1:32) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],9,32)-1)+1i*(2*randi([0,1],9,32)-1));
        X_DD_u4(26:32,9:29) = sqrt((1/2)*SNR(snr_iter))*((2*randi([0,1],7,21)-1)+1i*(2*randi([0,1],7,21)-1));
        
%         
        X_TF_u1 = F_m*X_DD_u1*F_n';
        X_TF_u2 = F_m*X_DD_u2*F_n';
        X_TF_u3 = F_m*X_DD_u3*F_n';
        X_TF_u4 = F_m*X_DD_u4*F_n';
        S_mat_u1 = Ptx*F_m'*X_TF_u1;
        S_mat_u2 = Ptx*F_m'*X_TF_u2;
        S_mat_u3 = Ptx*F_m'*X_TF_u3;
        S_mat_u4 = Ptx*F_m'*X_TF_u4;
        TxSamples_u1 = S_mat_u1(:).';
        TxSamples_u2 = S_mat_u2(:).';
        TxSamples_u3 = S_mat_u3(:).';
        TxSamples_u4 = S_mat_u4(:).';
        TxSamplesCP_u1 = [TxSamples_u1(M*N-Ncp_1+1:M*N) TxSamples_u1];
        TxSamplesCP_u2 = [TxSamples_u2(M*N-Ncp_2+1:M*N) TxSamples_u2];
        TxSamplesCP_u3 = [TxSamples_u3(M*N-Ncp_3+1:M*N) TxSamples_u3];
        TxSamplesCP_u4 = [TxSamples_u4(M*N-Ncp_4+1:M*N) TxSamples_u4];
        RxsamplesCP_u1 = 0;
        RxsamplesCP_u2 = 0;
        RxsamplesCP_u3 = 0;
        RxsamplesCP_u4 = 0;
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
        Rxsamples = RxsamplesCP_u1(Ncp_1+1:M*N+Ncp_1) + RxsamplesCP_u2(Ncp_2+1:M*N+Ncp_2) + RxsamplesCP_u3(Ncp_3+1:M*N+Ncp_3) + RxsamplesCP_u4(Ncp_4+1:M*N+Ncp_4) + noise;
        R_mat = reshape(Rxsamples.',M,N);
        Y_TF = F_m*Prx*R_mat;
        Y_DD = F_m'*Y_TF*F_n;
        
        % DD-domain CSI estimation
        [h_est_u1,delay_taps_est_u1,doppler_taps_est_u1,num_est_taps_u1,estimated_channel_matrix_u1] = Pilot_impulse_CE_modified(M,N,k_p_u1,l_p_u1,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u1 = H_eff_creation(M,N,num_est_taps_u1,h_est_u1,delay_taps_est_u1,doppler_taps_est_u1,Prx,Ptx);
        
        
        [h_est_u2,delay_taps_est_u2,doppler_taps_est_u2,num_est_taps_u2,estimated_channel_matrix_u2] = Pilot_impulse_CE_modified(M,N,k_p_u2,l_p_u2,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u2 = H_eff_creation(M,N,num_est_taps_u2,h_est_u2,delay_taps_est_u2,doppler_taps_est_u2,Prx,Ptx);

        [h_est_u3,delay_taps_est_u3,doppler_taps_est_u3,num_est_taps_u3,estimated_channel_matrix_u3] = Pilot_impulse_CE_modified(M,N,k_p_u3,l_p_u3,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u3 = H_eff_creation(M,N,num_est_taps_u3,h_est_u3,delay_taps_est_u3,doppler_taps_est_u3,Prx,Ptx);

        [h_est_u4,delay_taps_est_u4,doppler_taps_est_u4,num_est_taps_u4,estimated_channel_matrix_u4] = Pilot_impulse_CE_modified(M,N,k_p_u4,l_p_u4,l_max,k_max,thrld,...
                        Y_DD,pp_dB);
        H_eff_est_u4 = H_eff_creation(M,N,num_est_taps_u4,h_est_u4,delay_taps_est_u4,doppler_taps_est_u4,Prx,Ptx);

        
        %actual channel estimation
        % mse_impulse_pilot_u1(1,snr_iter) = mse_impulse_pilot_u1(1,snr_iter)+(((norm(estimated_channel_matrix_u1-true_channel_matrix_1,'fro'))^2)/(norm(true_channel_matrix_1,'fro')^2));
        % mse_impulse_pilot_u2(1,snr_iter) = mse_impulse_pilot_u2(1,snr_iter)+(((norm(estimated_channel_matrix_u2-true_channel_matrix_2,'fro'))^2)/(norm(true_channel_matrix_2,'fro')^2));
        % mse_impulse_pilot_u3(1,snr_iter) = mse_impulse_pilot_u3(1,snr_iter)+(((norm(estimated_channel_matrix_u3-true_channel_matrix_3,'fro'))^2)/(norm(true_channel_matrix_3,'fro')^2));
        % mse_impulse_pilot_u4(1,snr_iter) = mse_impulse_pilot_u4(1,snr_iter)+(((norm(estimated_channel_matrix_u4-true_channel_matrix_4,'fro'))^2)/(norm(true_channel_matrix_4,'fro')^2));

        %Embedded-data Detection
          x_DD_u1 = X_DD_u1(:);
          x_DD_u2 = X_DD_u2(:);
          x_DD_u3 = X_DD_u3(:);
          x_DD_u4 = X_DD_u4(:);
          x_DD_comp = [x_DD_u1;x_DD_u2;x_DD_u3;x_DD_u4];
          set_S_comp = find(x_DD_comp);
          ind_u1 = find(set_S_comp == 1);
          set_S_comp(ind_u1) = [];
          ind_u2 = find(set_S_comp == 1033);
          set_S_comp(ind_u2) = [];
          ind_u3 = find(set_S_comp == 2177);
          set_S_comp(ind_u3) = [];
          ind_u4 = find(set_S_comp == 3209);
          set_S_comp(ind_u4) = [];
          x_DD_S_comp = x_DD_comp(set_S_comp);    %only data
          y_DD_comp = Y_DD(:); 

        %Using estimated CSI
        y_DD_S_comp = y_DD_comp - (sqrt(10^(pp_dB/10))*H_eff_est_u1(:,1)) - (sqrt(10^(pp_dB/10))*H_eff_est_u2(:,9))- (sqrt(10^(pp_dB/10))*H_eff_est_u3(:,129))- (sqrt(10^(pp_dB/10))*H_eff_est_u4(:,137));
        H_eff_est_comp = [H_eff_est_u1 H_eff_est_u2 H_eff_est_u3 H_eff_est_u4];
        H_eff_est_S_comp = H_eff_est_comp(:,set_S_comp);
        x_DD_S_hat_LMMSE_comp = inv(H_eff_est_S_comp'*H_eff_est_S_comp + eye(length(set_S_comp))/SNR(snr_iter))*H_eff_est_S_comp'*y_DD_S_comp;
        DecodedBits_LMMSE_real = (real(x_DD_S_hat_LMMSE_comp)>=0);
        DecodedBits_LMMSE_imag = (imag(x_DD_S_hat_LMMSE_comp)>=0);
        DecodedSyms_LMMSE_real = sqrt((1/2)*SNR(snr_iter))*(2*DecodedBits_LMMSE_real-1);
        DecodedSyms_LMMSE_imag = sqrt((1/2)*SNR(snr_iter))*(2*DecodedBits_LMMSE_imag-1);
        BER_MMSE_ICSI_LMMSE(snr_iter) = BER_MMSE_ICSI_LMMSE(snr_iter) + sum(DecodedSyms_LMMSE_real ~= real(x_DD_S_comp))+ sum(DecodedSyms_LMMSE_imag ~= imag(x_DD_S_comp));
       
        % Using perfect CSI
        y_DD_S_comp = y_DD_comp - H_eff_u1(:,1) - (sqrt(10^(pp_dB/10))*H_eff_u2(:,9))- (sqrt(10^(pp_dB/10))*H_eff_u3(:,129))- (pp(sqrt^(10_dB/10))*H_eff_u4(:,137));   %note that 496 is the case for 32*32 and pilot placed at 16... automate this for different M and N with different pilot positions...
        H_eff_comp = [H_eff_u1 H_eff_u2 H_eff_u3 H_eff_u4];       
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

% mse_impulse_pilot_u1 = mse_impulse_pilot_u1/iter;
% mse_impulse_pilot_u2 = mse_impulse_pilot_u2/iter;
% mse_impulse_pilot_u3 = mse_impulse_pilot_u3/iter;
% mse_impulse_pilot_u4 = mse_impulse_pilot_u4/iter;

% semilogy(SNRdB,BER_MMSE_ICSI_LMMSE_u1,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% hold on; 
% grid on; 
% axis tight;
% semilogy(SNRdB,BER_MMSE_PCSI_LMMSE_u1,'r-.o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
% title('BER v/s SNR with LMMSE');
% xlabel('SNR(dB)'); ylabel('BER');
% legend('Imperfect CSI','Perfect CSI')
% 
% figure;
% semilogy(SNRdB,BER_MMSE_ICSI_ZF_u1,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% hold on; 
% grid on; 
% axis tight;
% semilogy(SNRdB,BER_MMSE_PCSI_ZF_u1,'r-.o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
% title('BER v/s SNR with ZF');
% xlabel('SNR(dB)'); ylabel('BER');
% legend('Imperfect CSI','Perfect CSI')
% 
% figure;
% semilogy(SNRdB,BER_MMSE_ICSI_LMMSE_u2,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% hold on; 
% grid on; 
% axis tight;
% semilogy(SNRdB,BER_MMSE_PCSI_LMMSE_u2,'r-.o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
% title('BER v/s SNR with LMMSE');
% xlabel('SNR(dB)'); ylabel('BER');
% legend('Imperfect CSI','Perfect CSI')
% 
% figure;
% semilogy(SNRdB,BER_MMSE_ICSI_ZF_u2,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% hold on; 
% grid on; 
% axis tight;
% semilogy(SNRdB,BER_MMSE_PCSI_ZF_u2,'r-.o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
% title('BER v/s SNR with ZF');
% xlabel('SNR(dB)'); ylabel('BER');
% legend('Imperfect CSI','Perfect CSI')
% 
% figure;
% semilogy(SNRdB,mse_impulse_pilot_u1,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% grid on; 
% axis tight;
% title('NMSE v/s SNR for 25dB pilot power');
% xlabel('SNR(dB)'); ylabel('NMSE');
% legend('Impulse Pilot Channel Estimation')
% 
% figure;
% semilogy(SNRdB,mse_impulse_pilot_u2,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% grid on; 
% axis tight;
% title('NMSE v/s SNR for 25dB pilot power');
% xlabel('SNR(dB)'); ylabel('NMSE');
% legend('Impulse Pilot Channel Estimation')
% 
% figure;
% semilogy(SNRdB,mse_impulse_pilot_u3,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% grid on; 
% axis tight;
% title('NMSE v/s SNR for 25dB pilot power');
% xlabel('SNR(dB)'); ylabel('NMSE');
% legend('Impulse Pilot Channel Estimation')
% 
% figure;
% semilogy(SNRdB,mse_impulse_pilot_u4,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
% grid on; 
% axis tight;
% title('NMSE v/s SNR for 25dB pilot power');
% xlabel('SNR(dB)'); ylabel('NMSE');
% legend('Impulse Pilot Channel Estimation')

figure;
semilogy(SNRdB,BER_MMSE_ICSI_LMMSE,'b-s','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
hold on;
semilogy(SNRdB,BER_MMSE_PCSI_LMMSE,'b-.o','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
hold on;
grid on;
axis tight;
%title('LMMSE VS ZF');
xlabel('SNR(dB)');
ylabel('BER');
% legend('Imperfect LMMSE U1','Perfect LMMSE U1','Imperfect LMMSE U2','Perfect LMMSE U2','Imperfect ZF U1','Perfect ZF U1','Imperfect ZF U2','Perfect ZF U2');
% legend('Imperfect LMMSE U1','Perfect LMMSE U1','Imperfect ZF U1','Perfect ZF U1');
legend('Imperfect CSI-LMMSE','Perfect CSI-LMMSE')
