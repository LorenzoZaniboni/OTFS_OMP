
close all;
clc;
clear all;

SNR_dB = -10:5:40;
SNR = 10.^(SNR_dB/10);
M = 32;  %grid size
N = 32;  %grid size
k_max = 4;   % max doppler tap
l_max = 8;   % max delay tap
G_tx = eye(M);  %transmit pusle shaping
G_rx = eye(N);  %receiver pulse shaping
taps = 4;       %no. of multipaths
iter = 1000;
F_m = 1/sqrt(M)*dftmtx(M);      %DFT matrix
F_n = 1/sqrt(N)*dftmtx(N);
X_dd_user_1 = zeros(M,N);    %grid of user 1 in DD domain
X_dd_user_2 = zeros(M,N);   %grid of user 2 in DD domain
X_dd_user_3 = zeros(M,N);   %grid of user 3 in DD domain
X_dd_user_4 = zeros(M,N);   %grid of user 4 in DD domain
X_dd_user_5 = zeros(M,N);   %grid of user 5 in DD domain
X_dd_user_6 = zeros(M,N);   %grid of user 6 in DD domain



l_p_u1 = 1;     %pilot position for user 1
k_p_u1 = 1;     %pilot position for user 1
l_p_u2 = 1;     %pilot position for user 2
k_p_u2 = 6;     %pilot position for user 2
l_p_u3 = 1;     %pilot position for user 3
k_p_u3 = 11;     %pilot position for user 3
l_p_u4 = 1;     %pilot position for user 4
k_p_u4 = 16;     %pilot position for user 4
l_p_u5 = 1;     %pilot position for user 5
k_p_u5 = 21;     %pilot position for user 5
l_p_u6 = 1;     %pilot position for user 6
k_p_u6 = 26;     %pilot position for user 6

X_dd_user_1(l_p_u1,k_p_u1) = 1;       % pilot impulse
X_dd_user_2(l_p_u2,k_p_u2) = 1;       % pilot impulse
X_dd_user_3(l_p_u3,k_p_u3) = 1;       % pilot impulse
X_dd_user_4(l_p_u4,k_p_u4) = 1;       % pilot impulse
X_dd_user_5(l_p_u5,k_p_u5) = 1;       % pilot impulse
X_dd_user_6(l_p_u6,k_p_u6) = 1;       % pilot impulse


% threshold = 1:2:7;
threshold = 3;
mse_impulse_pilot_user_1 = zeros(length(threshold),length(SNR));   %mean square error 
mse_impulse_pilot_user_2 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_3 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_4 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_5 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_6 = zeros(length(threshold),length(SNR));


mse_impulse_pilot_6_users = zeros(length(threshold),length(SNR));
    
    for n = 1:iter
    
    delay_taps_user_1 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_1 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_2 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_2 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_3 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_3 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_4 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_4 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_5 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_5 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_6 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_6 = randperm(k_max,taps)-1;   %to get zero as a possible tap



    % Cyclic prefix
    Ncp_user_1 = max(delay_taps_user_1);     %user 1 cyclic prefix
    Ncp_user_2 = max(delay_taps_user_2);    %user 2 cyclic prefix
    Ncp_user_3 = max(delay_taps_user_3);    %user 2 cyclic prefix
    Ncp_user_4 = max(delay_taps_user_4);    %user 2 cyclic prefix
    Ncp_user_5 = max(delay_taps_user_5);    %user 2 cyclic prefix
    Ncp_user_6 = max(delay_taps_user_6);    %user 2 cyclic prefix



    % Channel Coefficients
    h_user_1 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));  
    h_user_2 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    
    h_user_3 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    
    h_user_4 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    
    h_user_5 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    
    h_user_6 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    

    
    
    
    true_channel_matrix_user_1 = zeros(l_max,k_max);      %true matrix initialization
    true_channel_matrix_user_2 = zeros(l_max,k_max);
    true_channel_matrix_user_3 = zeros(l_max,k_max);
    true_channel_matrix_user_4 = zeros(l_max,k_max);
    true_channel_matrix_user_5 = zeros(l_max,k_max);
    true_channel_matrix_user_6 = zeros(l_max,k_max);




    for t = 1:taps    %creating the actual channel matrix in DD plane
        true_channel_matrix_user_1(delay_taps_user_1(t)+1,doppler_taps_user_1(t)+1) = h_user_1(t);  %since zero tap is considered, we need (0,0) possibility and is achieved by incrementing
        true_channel_matrix_user_2(delay_taps_user_2(t)+1,doppler_taps_user_2(t)+1) = h_user_2(t);
        true_channel_matrix_user_3(delay_taps_user_3(t)+1,doppler_taps_user_3(t)+1) = h_user_3(t);
        true_channel_matrix_user_4(delay_taps_user_4(t)+1,doppler_taps_user_4(t)+1) = h_user_4(t);
        true_channel_matrix_user_5(delay_taps_user_5(t)+1,doppler_taps_user_5(t)+1) = h_user_5(t);
        true_channel_matrix_user_6(delay_taps_user_6(t)+1,doppler_taps_user_6(t)+1) = h_user_6(t);

    end
    
    noise = sqrt(1/2)*(randn(1,M*N) + 1j*randn(1,M*N));

    for k = 1:length(SNR)
        disp("# of interation: " + n +  "; SNR: " + SNR_dB(k))
        X_DD_user_1 = sqrt(SNR(k))*X_dd_user_1;
        X_DD_user_2 = sqrt(SNR(k))*X_dd_user_2;
        X_DD_user_3 = sqrt(SNR(k))*X_dd_user_3;
        X_DD_user_4 = sqrt(SNR(k))*X_dd_user_4;
        X_DD_user_5 = sqrt(SNR(k))*X_dd_user_5;
        X_DD_user_6 = sqrt(SNR(k))*X_dd_user_6;


        X_TF_user_1 = F_m*X_DD_user_1*F_n';   %time frequency
        X_TF_user_2 = F_m*X_DD_user_2*F_n';
        X_TF_user_3 = F_m*X_DD_user_3*F_n';
        X_TF_user_4 = F_m*X_DD_user_4*F_n';
        X_TF_user_5 = F_m*X_DD_user_5*F_n';
        X_TF_user_6 = F_m*X_DD_user_6*F_n';

        S_user_1 = G_tx*F_m'*X_TF_user_1;       %heisenberg tranform
        S_user_2 = G_tx*F_m'*X_TF_user_2;
        S_user_3 = G_tx*F_m'*X_TF_user_3;
        S_user_4 = G_tx*F_m'*X_TF_user_4;
        S_user_5 = G_tx*F_m'*X_TF_user_5;
        S_user_6 = G_tx*F_m'*X_TF_user_6;

        
        TxSamples_user_1 = S_user_1(:).';       %time sequence
        TxSamples_user_2 = S_user_2(:).';
        TxSamples_user_3 = S_user_3(:).';
        TxSamples_user_4 = S_user_4(:).';
        TxSamples_user_5 = S_user_5(:).';
        TxSamples_user_6 = S_user_6(:).';

        TxSamplesCP_user_1 = [TxSamples_user_1(M*N-Ncp_user_1+1:M*N) TxSamples_user_1];  %adding CP
        TxSamplesCP_user_2 = [TxSamples_user_2(M*N-Ncp_user_2+1:M*N) TxSamples_user_2];
        TxSamplesCP_user_3 = [TxSamples_user_3(M*N-Ncp_user_3+1:M*N) TxSamples_user_3];
        TxSamplesCP_user_4 = [TxSamples_user_4(M*N-Ncp_user_4+1:M*N) TxSamples_user_4];
        TxSamplesCP_user_5 = [TxSamples_user_5(M*N-Ncp_user_5+1:M*N) TxSamples_user_5];
        TxSamplesCP_user_6 = [TxSamples_user_6(M*N-Ncp_user_6+1:M*N) TxSamples_user_6];

        
        RxsamplesCP_user_1 = zeros(1,M*N+Ncp_user_1)+1j*zeros(1,M*N+Ncp_user_1);   %initialzing received vector
        RxsamplesCP_user_2 = zeros(1,M*N+Ncp_user_2)+1j*zeros(1,M*N+Ncp_user_2);
        RxsamplesCP_user_3 = zeros(1,M*N+Ncp_user_3)+1j*zeros(1,M*N+Ncp_user_3);
        RxsamplesCP_user_4 = zeros(1,M*N+Ncp_user_4)+1j*zeros(1,M*N+Ncp_user_4);
        RxsamplesCP_user_5 = zeros(1,M*N+Ncp_user_5)+1j*zeros(1,M*N+Ncp_user_5);
        RxsamplesCP_user_6 = zeros(1,M*N+Ncp_user_6)+1j*zeros(1,M*N+Ncp_user_6);

        
        for tap_iter = 1:taps
            RxsamplesCP_user_1 = RxsamplesCP_user_1 + h_user_1(tap_iter)*circshift(TxSamplesCP_user_1.*exp(1j*2*pi/M*...
                (-Ncp_user_1:M*N-1)*doppler_taps_user_1(tap_iter)/N),[0,delay_taps_user_1(tap_iter)]);    %convolution
            RxsamplesCP_user_2 = RxsamplesCP_user_2 + h_user_2(tap_iter)*circshift(TxSamplesCP_user_2.*exp(1j*2*pi/M*...
                (-Ncp_user_2:M*N-1)*doppler_taps_user_2(tap_iter)/N),[0,delay_taps_user_2(tap_iter)]);
            RxsamplesCP_user_3 = RxsamplesCP_user_3 + h_user_3(tap_iter)*circshift(TxSamplesCP_user_3.*exp(1j*2*pi/M*...
                (-Ncp_user_3:M*N-1)*doppler_taps_user_3(tap_iter)/N),[0,delay_taps_user_3(tap_iter)]);
            RxsamplesCP_user_4 = RxsamplesCP_user_4 + h_user_4(tap_iter)*circshift(TxSamplesCP_user_4.*exp(1j*2*pi/M*...
                (-Ncp_user_4:M*N-1)*doppler_taps_user_4(tap_iter)/N),[0,delay_taps_user_4(tap_iter)]);
        
            RxsamplesCP_user_5 = RxsamplesCP_user_5 + h_user_5(tap_iter)*circshift(TxSamplesCP_user_5.*exp(1j*2*pi/M*...
                (-Ncp_user_5:M*N-1)*doppler_taps_user_5(tap_iter)/N),[0,delay_taps_user_5(tap_iter)]);
        
            RxsamplesCP_user_6 = RxsamplesCP_user_6 + h_user_6(tap_iter)*circshift(TxSamplesCP_user_6.*exp(1j*2*pi/M*...
                (-Ncp_user_6:M*N-1)*doppler_taps_user_6(tap_iter)/N),[0,delay_taps_user_6(tap_iter)]);
        

        end
        Rxsamples_user_1 = RxsamplesCP_user_1(Ncp_user_1+1:M*N+Ncp_user_1);  %removing CP
        Rxsamples_user_2 = RxsamplesCP_user_2(Ncp_user_2+1:M*N+Ncp_user_2);     
        Rxsamples_user_3 = RxsamplesCP_user_3(Ncp_user_3+1:M*N+Ncp_user_3);     
        Rxsamples_user_4 = RxsamplesCP_user_4(Ncp_user_4+1:M*N+Ncp_user_4);     
        Rxsamples_user_5 = RxsamplesCP_user_5(Ncp_user_5+1:M*N+Ncp_user_5);     
        Rxsamples_user_6 = RxsamplesCP_user_6(Ncp_user_6+1:M*N+Ncp_user_6);     



        Rxsamples = Rxsamples_user_1 + Rxsamples_user_2 + Rxsamples_user_3 + Rxsamples_user_4 + Rxsamples_user_5 + Rxsamples_user_6 + noise;   %adding noise
        R = reshape(Rxsamples.',M,N);     %back to grid of M*N
        Y_TF = F_m*G_rx*R;           %back to Time frequency
        Y_DD = F_m'*Y_TF*F_n;   %the delay doppler domain of the received vector

        % delay doppler and channel parameter estimation:Threshold 
        for th = 1:length(threshold)

            estimated_channel_matrix_u1 = special_impulse_search_function(M,N,l_p_u1,k_p_u1,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u2 = special_impulse_search_function(M,N,l_p_u2,k_p_u2,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u3 = special_impulse_search_function(M,N,l_p_u3,k_p_u3,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u4 = special_impulse_search_function(M,N,l_p_u4,k_p_u4,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u5 = special_impulse_search_function(M,N,l_p_u5,k_p_u5,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u6 = special_impulse_search_function(M,N,l_p_u6,k_p_u6,l_max,k_max,threshold(th),SNR(k),Y_DD);




            estimated_channel_vector = [estimated_channel_matrix_u1(:);estimated_channel_matrix_u2(:);estimated_channel_matrix_u3(:);estimated_channel_matrix_u4(:);estimated_channel_matrix_u5(:);estimated_channel_matrix_u6(:)];
            true_channel_vector = [true_channel_matrix_user_1(:);true_channel_matrix_user_2(:);true_channel_matrix_user_3(:);true_channel_matrix_user_4(:);true_channel_matrix_user_5(:);true_channel_matrix_user_6(:)];


            %MSE computation
            mse_impulse_pilot_user_1(th,k) = mse_impulse_pilot_user_1(th,k)+(((norm(estimated_channel_matrix_u1-true_channel_matrix_user_1,'fro'))^2)/(norm(true_channel_matrix_user_1,'fro')^2));
            mse_impulse_pilot_user_2(th,k) = mse_impulse_pilot_user_2(th,k)+(((norm(estimated_channel_matrix_u2-true_channel_matrix_user_2,'fro'))^2)/(norm(true_channel_matrix_user_2,'fro')^2));
            mse_impulse_pilot_user_3(th,k) = mse_impulse_pilot_user_3(th,k)+(((norm(estimated_channel_matrix_u3-true_channel_matrix_user_3,'fro'))^2)/(norm(true_channel_matrix_user_3,'fro')^2));
            mse_impulse_pilot_user_4(th,k) = mse_impulse_pilot_user_4(th,k)+(((norm(estimated_channel_matrix_u4-true_channel_matrix_user_4,'fro'))^2)/(norm(true_channel_matrix_user_4,'fro')^2));
            mse_impulse_pilot_user_5(th,k) = mse_impulse_pilot_user_5(th,k)+(((norm(estimated_channel_matrix_u5-true_channel_matrix_user_5,'fro'))^2)/(norm(true_channel_matrix_user_5,'fro')^2));
            mse_impulse_pilot_user_6(th,k) = mse_impulse_pilot_user_6(th,k)+(((norm(estimated_channel_matrix_u6-true_channel_matrix_user_6,'fro'))^2)/(norm(true_channel_matrix_user_6,'fro')^2));

             mse_impulse_pilot_6_users(th,k) = mse_impulse_pilot_6_users(th,k) + (((norm(estimated_channel_vector-true_channel_vector,'fro'))^2)/(norm(true_channel_vector,'fro')^2));  
        

        end
    end
    end



mse_impulse_pilot_6_users = mse_impulse_pilot_6_users/iter;



figure;
semilogy(SNR_dB,mse_impulse_pilot_6_users(1,:),'b-*','LineWidth',3,'MarkerFaceColor','b','MarkerSize',9.0);
hold on; grid on; 
axis tight;
title('6 Users Pilot-Impulse based CSI estimation with threshold 3');
xlabel('SNR(dB)'); ylabel('NMSE');
