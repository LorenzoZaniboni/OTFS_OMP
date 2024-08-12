
close all;
clc;
clear all;

SNR_dB = -10:5:40;
% SNR_dB = 40;
SNR = 10.^(SNR_dB/10);
% SNR = 100;
M = 32;  %grid size
N = 32;  %grid size
k_max = 5;   % max doppler tap
l_max = 7;   % max delay tap
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
X_dd_user_7 = zeros(M,N);   %grid of user 7 in DD domain
X_dd_user_8 = zeros(M,N);   %grid of user 8 in DD domain
X_dd_user_9 = zeros(M,N);   %grid of user 9 in DD domain
X_dd_user_10 = zeros(M,N);   %grid of user 10 in DD domain
X_dd_user_11 = zeros(M,N);    %grid of user 1 in DD domain
X_dd_user_12 = zeros(M,N);   %grid of user 2 in DD domain
X_dd_user_13 = zeros(M,N);   %grid of user 3 in DD domain
X_dd_user_14 = zeros(M,N);   %grid of user 4 in DD domain
X_dd_user_15 = zeros(M,N);   %grid of user 5 in DD domain
X_dd_user_16 = zeros(M,N);   %grid of user 6 in DD domain
X_dd_user_17 = zeros(M,N);   %grid of user 7 in DD domain
X_dd_user_18 = zeros(M,N);   %grid of user 8 in DD domain
X_dd_user_19 = zeros(M,N);   %grid of user 9 in DD domain
X_dd_user_20 = zeros(M,N);   %grid of user 10 in DD domain
X_dd_user_21 = zeros(M,N);   %grid of user 7 in DD domain
X_dd_user_22 = zeros(M,N);   %grid of user 8 in DD domain
X_dd_user_23 = zeros(M,N);   %grid of user 9 in DD domain
X_dd_user_24 = zeros(M,N);   %grid of user 10 in DD domain


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
l_p_u7 = 8;     %pilot position for user 7
k_p_u7 = 1;     %pilot position for user 7
l_p_u8 = 8;     %pilot position for user 8
k_p_u8 = 6;     %pilot position for user 8
l_p_u9 = 8;     %pilot position for user 9
k_p_u9 = 11;     %pilot position for user 9
l_p_u10 = 8;     %pilot position for user 10
k_p_u10 = 16;     %pilot position for user 10
l_p_u11 = 8;     %pilot position for user 1
k_p_u11 = 21;     %pilot position for user 1
l_p_u12 = 8;     %pilot position for user 2
k_p_u12 = 26;     %pilot position for user 2
l_p_u13 = 15;     %pilot position for user 3
k_p_u13 = 1;     %pilot position for user 3
l_p_u14 = 15;     %pilot position for user 4
k_p_u14 = 6;     %pilot position for user 4
l_p_u15 = 15;     %pilot position for user 5
k_p_u15 = 11;     %pilot position for user 5
l_p_u16 = 15;     %pilot position for user 6
k_p_u16 = 16;     %pilot position for user 6
l_p_u17 = 15;     %pilot position for user 7
k_p_u17 = 21;     %pilot position for user 7
l_p_u18 = 15;     %pilot position for user 8
k_p_u18 = 26;     %pilot position for user 8
l_p_u19 = 22;     %pilot position for user 9
k_p_u19 = 1;     %pilot position for user 9
l_p_u20 = 22;     %pilot position for user 10
k_p_u20 = 6;     %pilot position for user 10
l_p_u21 = 22;     %pilot position for user 7
k_p_u21 = 11;
l_p_u22 = 22;     %pilot position for user 7
k_p_u22 = 16;     %pilot position for user 8
l_p_u23 = 22;     %pilot position for user 8
k_p_u23 = 21;     %pilot position for user 9
l_p_u24 = 22;     %pilot position for user 9
k_p_u24 = 26;

X_dd_user_1(l_p_u1,k_p_u1) = 1;       % pilot impulse
X_dd_user_2(l_p_u2,k_p_u2) = 1;       % pilot impulse
X_dd_user_3(l_p_u3,k_p_u3) = 1;       % pilot impulse
X_dd_user_4(l_p_u4,k_p_u4) = 1;       % pilot impulse
X_dd_user_5(l_p_u5,k_p_u5) = 1;       % pilot impulse
X_dd_user_6(l_p_u6,k_p_u6) = 1;       % pilot impulse
X_dd_user_7(l_p_u7,k_p_u7) = 1;       % pilot impulse
X_dd_user_8(l_p_u8,k_p_u8) = 1;       % pilot impulse
X_dd_user_9(l_p_u9,k_p_u9) = 1;       % pilot impulse
X_dd_user_10(l_p_u10,k_p_u10) = 1;       % pilot impulse
X_dd_user_11(l_p_u11,k_p_u11) = 1;       % pilot impulse
X_dd_user_12(l_p_u12,k_p_u12) = 1;       % pilot impulse
X_dd_user_13(l_p_u13,k_p_u13) = 1;       % pilot impulse
X_dd_user_14(l_p_u14,k_p_u14) = 1;       % pilot impulse
X_dd_user_15(l_p_u15,k_p_u15) = 1;       % pilot impulse
X_dd_user_16(l_p_u16,k_p_u16) = 1;       % pilot impulse
X_dd_user_17(l_p_u17,k_p_u17) = 1;       % pilot impulse
X_dd_user_18(l_p_u18,k_p_u18) = 1;       % pilot impulse
X_dd_user_19(l_p_u19,k_p_u19) = 1;       % pilot impulse
X_dd_user_20(l_p_u20,k_p_u20) = 1;       % pilot impulse
X_dd_user_21(l_p_u21,k_p_u21) = 1;       % pilot impulse
X_dd_user_22(l_p_u22,k_p_u22) = 1;       % pilot impulse
X_dd_user_23(l_p_u23,k_p_u23) = 1;       % pilot impulse
X_dd_user_24(l_p_u24,k_p_u24) = 1;       % pilot impulse

% threshold = 1:2:7;
threshold = 3;
mse_impulse_pilot_user_1 = zeros(length(threshold),length(SNR));   %mean square error 
mse_impulse_pilot_user_2 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_3 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_4 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_5 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_6 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_7 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_8 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_9 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_10 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_11 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_12 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_13 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_14 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_15 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_16 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_17 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_18 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_19 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_20 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_21 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_22 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_23 = zeros(length(threshold),length(SNR));
mse_impulse_pilot_user_24 = zeros(length(threshold),length(SNR));
    
    for n = 1:iter
    n
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
    delay_taps_user_7 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_7 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_8 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_8 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_9 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_9 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_10 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_10 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_11 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_11 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_12 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_12 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_13 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_13 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_14 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_14 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_15 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_15 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_16 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_16 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_17 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_17 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_18 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_18 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_19 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_19 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_20 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_20 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_21 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_21 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_22 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_22 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_23 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_23 = randperm(k_max,taps)-1;   %to get zero as a possible tap
    delay_taps_user_24 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    doppler_taps_user_24 = randperm(k_max,taps)-1;   %to get zero as a possible tap



    Ncp_user_1 = max(delay_taps_user_1);     %user 1 cyclic prefix
    Ncp_user_2 = max(delay_taps_user_2);    %user 2 cyclic prefix
    Ncp_user_3 = max(delay_taps_user_3);    %user 2 cyclic prefix
    Ncp_user_4 = max(delay_taps_user_4);    %user 2 cyclic prefix
    Ncp_user_5 = max(delay_taps_user_5);    %user 2 cyclic prefix
    Ncp_user_6 = max(delay_taps_user_6);    %user 2 cyclic prefix
    Ncp_user_7 = max(delay_taps_user_7);    %user 2 cyclic prefix
    Ncp_user_8 = max(delay_taps_user_8);    %user 2 cyclic prefix
    Ncp_user_9 = max(delay_taps_user_9);    %user 2 cyclic prefix
    Ncp_user_10 = max(delay_taps_user_10);    %user 2 cyclic prefix
    Ncp_user_11 = max(delay_taps_user_11);     %user 1 cyclic prefix
    Ncp_user_12 = max(delay_taps_user_12);    %user 2 cyclic prefix
    Ncp_user_13 = max(delay_taps_user_13);    %user 2 cyclic prefix
    Ncp_user_14 = max(delay_taps_user_14);    %user 2 cyclic prefix
    Ncp_user_15 = max(delay_taps_user_15);    %user 2 cyclic prefix
    Ncp_user_16 = max(delay_taps_user_16);    %user 2 cyclic prefix
    Ncp_user_17 = max(delay_taps_user_17);    %user 2 cyclic prefix
    Ncp_user_18 = max(delay_taps_user_18);    %user 2 cyclic prefix
    Ncp_user_19 = max(delay_taps_user_19);    %user 2 cyclic prefix
    Ncp_user_20 = max(delay_taps_user_20);    %user 2 cyclic prefix
    Ncp_user_21 = max(delay_taps_user_21);    %user 2 cyclic prefix
    Ncp_user_22 = max(delay_taps_user_22);    %user 2 cyclic prefix
    Ncp_user_23 = max(delay_taps_user_23);    %user 2 cyclic prefix
    Ncp_user_24 = max(delay_taps_user_24);    %user 2 cyclic prefix



    h_user_1 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));   %channel coefficients of user 1
    h_user_2 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_3 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_4 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_5 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_6 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_7 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_8 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_9 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_10 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_11 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));   %channel coefficients of user 1
    h_user_12 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_13 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_14 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_15 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_16 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_17 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_18 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_19 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_20 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_21 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_22 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_23 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    h_user_24 = (1/sqrt(2))*(randn(1,taps) + 1i*randn(1,taps));    %channel coefficients of user 2
    
    
    
    true_channel_matrix_user_1 = zeros(l_max,k_max);      %true matrix initialization
    true_channel_matrix_user_2 = zeros(l_max,k_max);
    true_channel_matrix_user_3 = zeros(l_max,k_max);
    true_channel_matrix_user_4 = zeros(l_max,k_max);
    true_channel_matrix_user_5 = zeros(l_max,k_max);
    true_channel_matrix_user_6 = zeros(l_max,k_max);
    true_channel_matrix_user_7 = zeros(l_max,k_max);
    true_channel_matrix_user_8 = zeros(l_max,k_max);
    true_channel_matrix_user_9 = zeros(l_max,k_max);
    true_channel_matrix_user_10 = zeros(l_max,k_max);
    true_channel_matrix_user_11 = zeros(l_max,k_max);      %true matrix initialization
    true_channel_matrix_user_12 = zeros(l_max,k_max);
    true_channel_matrix_user_13 = zeros(l_max,k_max);
    true_channel_matrix_user_14 = zeros(l_max,k_max);
    true_channel_matrix_user_15 = zeros(l_max,k_max);
    true_channel_matrix_user_16 = zeros(l_max,k_max);
    true_channel_matrix_user_17 = zeros(l_max,k_max);
    true_channel_matrix_user_18 = zeros(l_max,k_max);
    true_channel_matrix_user_19 = zeros(l_max,k_max);
    true_channel_matrix_user_20 = zeros(l_max,k_max);
    true_channel_matrix_user_21 = zeros(l_max,k_max);
    true_channel_matrix_user_22 = zeros(l_max,k_max);
    true_channel_matrix_user_23 = zeros(l_max,k_max);
    true_channel_matrix_user_24 = zeros(l_max,k_max);



    for t = 1:taps    %creating the actual channel matrix in DD plane
        true_channel_matrix_user_1(delay_taps_user_1(t)+1,doppler_taps_user_1(t)+1) = h_user_1(t);  %since zero tap is considered, we need (0,0) possibility and is achieved by incrementing
        true_channel_matrix_user_2(delay_taps_user_2(t)+1,doppler_taps_user_2(t)+1) = h_user_2(t);
        true_channel_matrix_user_3(delay_taps_user_3(t)+1,doppler_taps_user_3(t)+1) = h_user_3(t);
        true_channel_matrix_user_4(delay_taps_user_4(t)+1,doppler_taps_user_4(t)+1) = h_user_4(t);
        true_channel_matrix_user_5(delay_taps_user_5(t)+1,doppler_taps_user_5(t)+1) = h_user_5(t);
        true_channel_matrix_user_6(delay_taps_user_6(t)+1,doppler_taps_user_6(t)+1) = h_user_6(t);
        true_channel_matrix_user_7(delay_taps_user_7(t)+1,doppler_taps_user_7(t)+1) = h_user_7(t);
        true_channel_matrix_user_8(delay_taps_user_8(t)+1,doppler_taps_user_8(t)+1) = h_user_8(t);
        true_channel_matrix_user_9(delay_taps_user_9(t)+1,doppler_taps_user_9(t)+1) = h_user_9(t);
        true_channel_matrix_user_10(delay_taps_user_10(t)+1,doppler_taps_user_10(t)+1) = h_user_10(t);
        true_channel_matrix_user_11(delay_taps_user_11(t)+1,doppler_taps_user_11(t)+1) = h_user_11(t);  %since zero tap is considered, we need (0,0) possibility and is achieved by incrementing
        true_channel_matrix_user_12(delay_taps_user_12(t)+1,doppler_taps_user_12(t)+1) = h_user_12(t);
        true_channel_matrix_user_13(delay_taps_user_13(t)+1,doppler_taps_user_13(t)+1) = h_user_13(t);
        true_channel_matrix_user_14(delay_taps_user_14(t)+1,doppler_taps_user_14(t)+1) = h_user_14(t);
        true_channel_matrix_user_15(delay_taps_user_15(t)+1,doppler_taps_user_15(t)+1) = h_user_15(t);
        true_channel_matrix_user_16(delay_taps_user_16(t)+1,doppler_taps_user_16(t)+1) = h_user_16(t);
        true_channel_matrix_user_17(delay_taps_user_17(t)+1,doppler_taps_user_17(t)+1) = h_user_17(t);
        true_channel_matrix_user_18(delay_taps_user_18(t)+1,doppler_taps_user_18(t)+1) = h_user_18(t);
        true_channel_matrix_user_19(delay_taps_user_19(t)+1,doppler_taps_user_19(t)+1) = h_user_19(t);
        true_channel_matrix_user_20(delay_taps_user_20(t)+1,doppler_taps_user_20(t)+1) = h_user_20(t);
        true_channel_matrix_user_21(delay_taps_user_21(t)+1,doppler_taps_user_21(t)+1) = h_user_21(t);
        true_channel_matrix_user_22(delay_taps_user_22(t)+1,doppler_taps_user_22(t)+1) = h_user_22(t);
        true_channel_matrix_user_23(delay_taps_user_23(t)+1,doppler_taps_user_23(t)+1) = h_user_23(t);
        true_channel_matrix_user_24(delay_taps_user_24(t)+1,doppler_taps_user_24(t)+1) = h_user_24(t);
    end
    
    noise = sqrt(1/2)*(randn(1,M*N) + 1j*randn(1,M*N));

    for k = 1:length(SNR)

        X_DD_user_1 = sqrt(SNR(k))*X_dd_user_1;
        X_DD_user_2 = sqrt(SNR(k))*X_dd_user_2;
        X_DD_user_3 = sqrt(SNR(k))*X_dd_user_3;
        X_DD_user_4 = sqrt(SNR(k))*X_dd_user_4;
        X_DD_user_5 = sqrt(SNR(k))*X_dd_user_5;
        X_DD_user_6 = sqrt(SNR(k))*X_dd_user_6;
        X_DD_user_7 = sqrt(SNR(k))*X_dd_user_7;
        X_DD_user_8 = sqrt(SNR(k))*X_dd_user_8;
        X_DD_user_9 = sqrt(SNR(k))*X_dd_user_9;
        X_DD_user_10 = sqrt(SNR(k))*X_dd_user_10;
        X_DD_user_11 = sqrt(SNR(k))*X_dd_user_11;
        X_DD_user_12 = sqrt(SNR(k))*X_dd_user_12;
        X_DD_user_13 = sqrt(SNR(k))*X_dd_user_13;
        X_DD_user_14 = sqrt(SNR(k))*X_dd_user_14;
        X_DD_user_15 = sqrt(SNR(k))*X_dd_user_15;
        X_DD_user_16 = sqrt(SNR(k))*X_dd_user_16;
        X_DD_user_17 = sqrt(SNR(k))*X_dd_user_17;
        X_DD_user_18 = sqrt(SNR(k))*X_dd_user_18;
        X_DD_user_19 = sqrt(SNR(k))*X_dd_user_19;
        X_DD_user_20 = sqrt(SNR(k))*X_dd_user_20;
        X_DD_user_21 = sqrt(SNR(k))*X_dd_user_21;
        X_DD_user_22 = sqrt(SNR(k))*X_dd_user_22;
        X_DD_user_23 = sqrt(SNR(k))*X_dd_user_23;
        X_DD_user_24 = sqrt(SNR(k))*X_dd_user_24;


        X_TF_user_1 = F_m*X_DD_user_1*F_n';   %time frequency
        X_TF_user_2 = F_m*X_DD_user_2*F_n';
        X_TF_user_3 = F_m*X_DD_user_3*F_n';
        X_TF_user_4 = F_m*X_DD_user_4*F_n';
        X_TF_user_5 = F_m*X_DD_user_5*F_n';
        X_TF_user_6 = F_m*X_DD_user_6*F_n';
        X_TF_user_7 = F_m*X_DD_user_7*F_n';
        X_TF_user_8 = F_m*X_DD_user_8*F_n';
        X_TF_user_9 = F_m*X_DD_user_9*F_n';
        X_TF_user_10 = F_m*X_DD_user_10*F_n';
        X_TF_user_11 = F_m*X_DD_user_11*F_n';   %time frequency
        X_TF_user_12 = F_m*X_DD_user_12*F_n';
        X_TF_user_13 = F_m*X_DD_user_13*F_n';
        X_TF_user_14 = F_m*X_DD_user_14*F_n';
        X_TF_user_15 = F_m*X_DD_user_15*F_n';
        X_TF_user_16 = F_m*X_DD_user_16*F_n';
        X_TF_user_17 = F_m*X_DD_user_17*F_n';
        X_TF_user_18 = F_m*X_DD_user_18*F_n';
        X_TF_user_19 = F_m*X_DD_user_19*F_n';
        X_TF_user_20 = F_m*X_DD_user_20*F_n';
        X_TF_user_21 = F_m*X_DD_user_21*F_n';
        X_TF_user_22 = F_m*X_DD_user_22*F_n';
        X_TF_user_23 = F_m*X_DD_user_23*F_n';
        X_TF_user_24 = F_m*X_DD_user_24*F_n';

        S_user_1 = G_tx*F_m'*X_TF_user_1;       %heisenberg tranform
        S_user_2 = G_tx*F_m'*X_TF_user_2;
        S_user_3 = G_tx*F_m'*X_TF_user_3;
        S_user_4 = G_tx*F_m'*X_TF_user_4;
        S_user_5 = G_tx*F_m'*X_TF_user_5;
        S_user_6 = G_tx*F_m'*X_TF_user_6;
        S_user_7 = G_tx*F_m'*X_TF_user_7;
        S_user_8 = G_tx*F_m'*X_TF_user_8;
        S_user_9 = G_tx*F_m'*X_TF_user_9;
        S_user_10 = G_tx*F_m'*X_TF_user_10;
        S_user_11 = G_tx*F_m'*X_TF_user_11;       %heisenberg tranform
        S_user_12 = G_tx*F_m'*X_TF_user_12;
        S_user_13 = G_tx*F_m'*X_TF_user_13;
        S_user_14 = G_tx*F_m'*X_TF_user_14;
        S_user_15 = G_tx*F_m'*X_TF_user_15;
        S_user_16 = G_tx*F_m'*X_TF_user_16;
        S_user_17 = G_tx*F_m'*X_TF_user_17;
        S_user_18 = G_tx*F_m'*X_TF_user_18;
        S_user_19 = G_tx*F_m'*X_TF_user_19;
        S_user_20 = G_tx*F_m'*X_TF_user_20;
        S_user_21 = G_tx*F_m'*X_TF_user_21;
        S_user_22 = G_tx*F_m'*X_TF_user_22;
        S_user_23 = G_tx*F_m'*X_TF_user_23;
        S_user_24 = G_tx*F_m'*X_TF_user_24;

        
        TxSamples_user_1 = S_user_1(:).';       %time sequence
        TxSamples_user_2 = S_user_2(:).';
        TxSamples_user_3 = S_user_3(:).';
        TxSamples_user_4 = S_user_4(:).';
        TxSamples_user_5 = S_user_5(:).';
        TxSamples_user_6 = S_user_6(:).';
        TxSamples_user_7 = S_user_7(:).';
        TxSamples_user_8 = S_user_8(:).';
        TxSamples_user_9 = S_user_9(:).';
        TxSamples_user_10 = S_user_10(:).';
        TxSamples_user_11 = S_user_11(:).';       %time sequence
        TxSamples_user_12 = S_user_12(:).';
        TxSamples_user_13 = S_user_13(:).';
        TxSamples_user_14 = S_user_14(:).';
        TxSamples_user_15 = S_user_15(:).';
        TxSamples_user_16 = S_user_16(:).';
        TxSamples_user_17 = S_user_17(:).';
        TxSamples_user_18 = S_user_18(:).';
        TxSamples_user_19 = S_user_19(:).';
        TxSamples_user_20 = S_user_20(:).';
        TxSamples_user_21 = S_user_21(:).';
        TxSamples_user_22 = S_user_22(:).';
        TxSamples_user_23 = S_user_23(:).';
        TxSamples_user_24 = S_user_24(:).';

        TxSamplesCP_user_1 = [TxSamples_user_1(M*N-Ncp_user_1+1:M*N) TxSamples_user_1];  %adding CP
        TxSamplesCP_user_2 = [TxSamples_user_2(M*N-Ncp_user_2+1:M*N) TxSamples_user_2];
        TxSamplesCP_user_3 = [TxSamples_user_3(M*N-Ncp_user_3+1:M*N) TxSamples_user_3];
        TxSamplesCP_user_4 = [TxSamples_user_4(M*N-Ncp_user_4+1:M*N) TxSamples_user_4];
        TxSamplesCP_user_5 = [TxSamples_user_5(M*N-Ncp_user_5+1:M*N) TxSamples_user_5];
        TxSamplesCP_user_6 = [TxSamples_user_6(M*N-Ncp_user_6+1:M*N) TxSamples_user_6];
        TxSamplesCP_user_7 = [TxSamples_user_7(M*N-Ncp_user_7+1:M*N) TxSamples_user_7];
        TxSamplesCP_user_8 = [TxSamples_user_8(M*N-Ncp_user_8+1:M*N) TxSamples_user_8];
        TxSamplesCP_user_9 = [TxSamples_user_9(M*N-Ncp_user_9+1:M*N) TxSamples_user_9];
        TxSamplesCP_user_10 = [TxSamples_user_10(M*N-Ncp_user_10+1:M*N) TxSamples_user_10];
        TxSamplesCP_user_11 = [TxSamples_user_11(M*N-Ncp_user_11+1:M*N) TxSamples_user_11];  %adding CP
        TxSamplesCP_user_12 = [TxSamples_user_12(M*N-Ncp_user_12+1:M*N) TxSamples_user_12];
        TxSamplesCP_user_13 = [TxSamples_user_13(M*N-Ncp_user_13+1:M*N) TxSamples_user_13];
        TxSamplesCP_user_14 = [TxSamples_user_14(M*N-Ncp_user_14+1:M*N) TxSamples_user_14];
        TxSamplesCP_user_15 = [TxSamples_user_15(M*N-Ncp_user_15+1:M*N) TxSamples_user_15];
        TxSamplesCP_user_16 = [TxSamples_user_16(M*N-Ncp_user_16+1:M*N) TxSamples_user_16];
        TxSamplesCP_user_17 = [TxSamples_user_17(M*N-Ncp_user_17+1:M*N) TxSamples_user_17];
        TxSamplesCP_user_18 = [TxSamples_user_18(M*N-Ncp_user_18+1:M*N) TxSamples_user_18];
        TxSamplesCP_user_19 = [TxSamples_user_19(M*N-Ncp_user_19+1:M*N) TxSamples_user_19];
        TxSamplesCP_user_20 = [TxSamples_user_20(M*N-Ncp_user_20+1:M*N) TxSamples_user_20];
        TxSamplesCP_user_21 = [TxSamples_user_21(M*N-Ncp_user_21+1:M*N) TxSamples_user_21];
        TxSamplesCP_user_22 = [TxSamples_user_22(M*N-Ncp_user_22+1:M*N) TxSamples_user_22];
        TxSamplesCP_user_23 = [TxSamples_user_23(M*N-Ncp_user_23+1:M*N) TxSamples_user_23];
        TxSamplesCP_user_24 = [TxSamples_user_24(M*N-Ncp_user_24+1:M*N) TxSamples_user_24];
        
        RxsamplesCP_user_1 = zeros(1,M*N+Ncp_user_1)+1j*zeros(1,M*N+Ncp_user_1);   %initialzing received vector
        RxsamplesCP_user_2 = zeros(1,M*N+Ncp_user_2)+1j*zeros(1,M*N+Ncp_user_2);
        RxsamplesCP_user_3 = zeros(1,M*N+Ncp_user_3)+1j*zeros(1,M*N+Ncp_user_3);
        RxsamplesCP_user_4 = zeros(1,M*N+Ncp_user_4)+1j*zeros(1,M*N+Ncp_user_4);
        RxsamplesCP_user_5 = zeros(1,M*N+Ncp_user_5)+1j*zeros(1,M*N+Ncp_user_5);
        RxsamplesCP_user_6 = zeros(1,M*N+Ncp_user_6)+1j*zeros(1,M*N+Ncp_user_6);
        RxsamplesCP_user_7 = zeros(1,M*N+Ncp_user_7)+1j*zeros(1,M*N+Ncp_user_7);
        RxsamplesCP_user_8 = zeros(1,M*N+Ncp_user_8)+1j*zeros(1,M*N+Ncp_user_8);
        RxsamplesCP_user_9 = zeros(1,M*N+Ncp_user_9)+1j*zeros(1,M*N+Ncp_user_9);
        RxsamplesCP_user_10 = zeros(1,M*N+Ncp_user_10)+1j*zeros(1,M*N+Ncp_user_10);
        RxsamplesCP_user_11 = zeros(1,M*N+Ncp_user_11)+1j*zeros(1,M*N+Ncp_user_11);   %initialzing received vector
        RxsamplesCP_user_12 = zeros(1,M*N+Ncp_user_12)+1j*zeros(1,M*N+Ncp_user_12);
        RxsamplesCP_user_13 = zeros(1,M*N+Ncp_user_13)+1j*zeros(1,M*N+Ncp_user_13);
        RxsamplesCP_user_14 = zeros(1,M*N+Ncp_user_14)+1j*zeros(1,M*N+Ncp_user_14);
        RxsamplesCP_user_15 = zeros(1,M*N+Ncp_user_15)+1j*zeros(1,M*N+Ncp_user_15);
        RxsamplesCP_user_16 = zeros(1,M*N+Ncp_user_16)+1j*zeros(1,M*N+Ncp_user_16);
        RxsamplesCP_user_17 = zeros(1,M*N+Ncp_user_17)+1j*zeros(1,M*N+Ncp_user_17);
        RxsamplesCP_user_18 = zeros(1,M*N+Ncp_user_18)+1j*zeros(1,M*N+Ncp_user_18);
        RxsamplesCP_user_19 = zeros(1,M*N+Ncp_user_19)+1j*zeros(1,M*N+Ncp_user_19);
        RxsamplesCP_user_20 = zeros(1,M*N+Ncp_user_20)+1j*zeros(1,M*N+Ncp_user_20);
        RxsamplesCP_user_21 = zeros(1,M*N+Ncp_user_21)+1j*zeros(1,M*N+Ncp_user_21);
        RxsamplesCP_user_22 = zeros(1,M*N+Ncp_user_22)+1j*zeros(1,M*N+Ncp_user_22);
        RxsamplesCP_user_23 = zeros(1,M*N+Ncp_user_23)+1j*zeros(1,M*N+Ncp_user_23);
        RxsamplesCP_user_24 = zeros(1,M*N+Ncp_user_24)+1j*zeros(1,M*N+Ncp_user_24);
        
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
        
            RxsamplesCP_user_7 = RxsamplesCP_user_7 + h_user_7(tap_iter)*circshift(TxSamplesCP_user_7.*exp(1j*2*pi/M*...
                (-Ncp_user_7:M*N-1)*doppler_taps_user_7(tap_iter)/N),[0,delay_taps_user_7(tap_iter)]);
        
            RxsamplesCP_user_8 = RxsamplesCP_user_8 + h_user_8(tap_iter)*circshift(TxSamplesCP_user_8.*exp(1j*2*pi/M*...
                (-Ncp_user_8:M*N-1)*doppler_taps_user_8(tap_iter)/N),[0,delay_taps_user_8(tap_iter)]);
        
            RxsamplesCP_user_9 = RxsamplesCP_user_9 + h_user_9(tap_iter)*circshift(TxSamplesCP_user_9.*exp(1j*2*pi/M*...
                (-Ncp_user_9:M*N-1)*doppler_taps_user_9(tap_iter)/N),[0,delay_taps_user_9(tap_iter)]);
        
            RxsamplesCP_user_10 = RxsamplesCP_user_10 + h_user_10(tap_iter)*circshift(TxSamplesCP_user_10.*exp(1j*2*pi/M*...
                (-Ncp_user_10:M*N-1)*doppler_taps_user_10(tap_iter)/N),[0,delay_taps_user_10(tap_iter)]);
            RxsamplesCP_user_11 = RxsamplesCP_user_11 + h_user_11(tap_iter)*circshift(TxSamplesCP_user_11.*exp(1j*2*pi/M*...
                (-Ncp_user_11:M*N-1)*doppler_taps_user_11(tap_iter)/N),[0,delay_taps_user_11(tap_iter)]);    %convolution
            RxsamplesCP_user_12 = RxsamplesCP_user_12 + h_user_12(tap_iter)*circshift(TxSamplesCP_user_12.*exp(1j*2*pi/M*...
                (-Ncp_user_12:M*N-1)*doppler_taps_user_12(tap_iter)/N),[0,delay_taps_user_12(tap_iter)]);
            RxsamplesCP_user_13 = RxsamplesCP_user_13 + h_user_13(tap_iter)*circshift(TxSamplesCP_user_13.*exp(1j*2*pi/M*...
                (-Ncp_user_13:M*N-1)*doppler_taps_user_13(tap_iter)/N),[0,delay_taps_user_13(tap_iter)]);
            RxsamplesCP_user_14 = RxsamplesCP_user_14 + h_user_14(tap_iter)*circshift(TxSamplesCP_user_14.*exp(1j*2*pi/M*...
                (-Ncp_user_14:M*N-1)*doppler_taps_user_14(tap_iter)/N),[0,delay_taps_user_14(tap_iter)]);
        
            RxsamplesCP_user_15 = RxsamplesCP_user_15 + h_user_15(tap_iter)*circshift(TxSamplesCP_user_15.*exp(1j*2*pi/M*...
                (-Ncp_user_15:M*N-1)*doppler_taps_user_15(tap_iter)/N),[0,delay_taps_user_15(tap_iter)]);
        
            RxsamplesCP_user_16 = RxsamplesCP_user_16 + h_user_16(tap_iter)*circshift(TxSamplesCP_user_16.*exp(1j*2*pi/M*...
                (-Ncp_user_16:M*N-1)*doppler_taps_user_16(tap_iter)/N),[0,delay_taps_user_16(tap_iter)]);
        
            RxsamplesCP_user_17 = RxsamplesCP_user_17 + h_user_17(tap_iter)*circshift(TxSamplesCP_user_17.*exp(1j*2*pi/M*...
                (-Ncp_user_17:M*N-1)*doppler_taps_user_17(tap_iter)/N),[0,delay_taps_user_17(tap_iter)]);
        
            RxsamplesCP_user_18 = RxsamplesCP_user_18 + h_user_18(tap_iter)*circshift(TxSamplesCP_user_18.*exp(1j*2*pi/M*...
                (-Ncp_user_18:M*N-1)*doppler_taps_user_18(tap_iter)/N),[0,delay_taps_user_18(tap_iter)]);
        
            RxsamplesCP_user_19 = RxsamplesCP_user_19 + h_user_19(tap_iter)*circshift(TxSamplesCP_user_19.*exp(1j*2*pi/M*...
                (-Ncp_user_19:M*N-1)*doppler_taps_user_19(tap_iter)/N),[0,delay_taps_user_19(tap_iter)]);
        
            RxsamplesCP_user_20 = RxsamplesCP_user_20 + h_user_20(tap_iter)*circshift(TxSamplesCP_user_20.*exp(1j*2*pi/M*...
                (-Ncp_user_20:M*N-1)*doppler_taps_user_20(tap_iter)/N),[0,delay_taps_user_20(tap_iter)]);
            RxsamplesCP_user_21 = RxsamplesCP_user_21 + h_user_21(tap_iter)*circshift(TxSamplesCP_user_21.*exp(1j*2*pi/M*...
                (-Ncp_user_21:M*N-1)*doppler_taps_user_21(tap_iter)/N),[0,delay_taps_user_21(tap_iter)]);
        
            RxsamplesCP_user_22 = RxsamplesCP_user_22 + h_user_22(tap_iter)*circshift(TxSamplesCP_user_22.*exp(1j*2*pi/M*...
                (-Ncp_user_22:M*N-1)*doppler_taps_user_22(tap_iter)/N),[0,delay_taps_user_22(tap_iter)]);
        
            RxsamplesCP_user_23 = RxsamplesCP_user_23 + h_user_23(tap_iter)*circshift(TxSamplesCP_user_23.*exp(1j*2*pi/M*...
                (-Ncp_user_23:M*N-1)*doppler_taps_user_23(tap_iter)/N),[0,delay_taps_user_23(tap_iter)]);
        
            RxsamplesCP_user_24 = RxsamplesCP_user_24 + h_user_24(tap_iter)*circshift(TxSamplesCP_user_24.*exp(1j*2*pi/M*...
                (-Ncp_user_24:M*N-1)*doppler_taps_user_24(tap_iter)/N),[0,delay_taps_user_24(tap_iter)]);

        end
        Rxsamples_user_1 = RxsamplesCP_user_1(Ncp_user_1+1:M*N+Ncp_user_1);  %removing CP
        Rxsamples_user_2 = RxsamplesCP_user_2(Ncp_user_2+1:M*N+Ncp_user_2);     
        Rxsamples_user_3 = RxsamplesCP_user_3(Ncp_user_3+1:M*N+Ncp_user_3);     
        Rxsamples_user_4 = RxsamplesCP_user_4(Ncp_user_4+1:M*N+Ncp_user_4);     
        Rxsamples_user_5 = RxsamplesCP_user_5(Ncp_user_5+1:M*N+Ncp_user_5);     
        Rxsamples_user_6 = RxsamplesCP_user_6(Ncp_user_6+1:M*N+Ncp_user_6);     
        Rxsamples_user_7 = RxsamplesCP_user_7(Ncp_user_7+1:M*N+Ncp_user_7);     
        Rxsamples_user_8 = RxsamplesCP_user_8(Ncp_user_8+1:M*N+Ncp_user_8);     
        Rxsamples_user_9 = RxsamplesCP_user_9(Ncp_user_9+1:M*N+Ncp_user_9);     
        Rxsamples_user_10 = RxsamplesCP_user_10(Ncp_user_10+1:M*N+Ncp_user_10);  
        Rxsamples_user_11 = RxsamplesCP_user_11(Ncp_user_11+1:M*N+Ncp_user_11);  %removing CP
        Rxsamples_user_12 = RxsamplesCP_user_12(Ncp_user_12+1:M*N+Ncp_user_12);     
        Rxsamples_user_13 = RxsamplesCP_user_13(Ncp_user_13+1:M*N+Ncp_user_13);     
        Rxsamples_user_14 = RxsamplesCP_user_14(Ncp_user_14+1:M*N+Ncp_user_14);     
        Rxsamples_user_15 = RxsamplesCP_user_15(Ncp_user_15+1:M*N+Ncp_user_15);     
        Rxsamples_user_16 = RxsamplesCP_user_16(Ncp_user_16+1:M*N+Ncp_user_16);     
        Rxsamples_user_17 = RxsamplesCP_user_17(Ncp_user_17+1:M*N+Ncp_user_17);     
        Rxsamples_user_18 = RxsamplesCP_user_18(Ncp_user_18+1:M*N+Ncp_user_18);     
        Rxsamples_user_19 = RxsamplesCP_user_19(Ncp_user_19+1:M*N+Ncp_user_19);     
        Rxsamples_user_20 = RxsamplesCP_user_20(Ncp_user_20+1:M*N+Ncp_user_20);  
        Rxsamples_user_21 = RxsamplesCP_user_21(Ncp_user_21+1:M*N+Ncp_user_21);     
        Rxsamples_user_22 = RxsamplesCP_user_22(Ncp_user_22+1:M*N+Ncp_user_22);     
        Rxsamples_user_23 = RxsamplesCP_user_23(Ncp_user_23+1:M*N+Ncp_user_23);     
        Rxsamples_user_24 = RxsamplesCP_user_24(Ncp_user_24+1:M*N+Ncp_user_24);  


        Rxsamples = Rxsamples_user_1 + Rxsamples_user_2 + Rxsamples_user_3 + Rxsamples_user_4 + Rxsamples_user_5 + Rxsamples_user_6 + Rxsamples_user_7 + Rxsamples_user_8 + Rxsamples_user_9 + Rxsamples_user_10 + Rxsamples_user_11 + Rxsamples_user_12 + Rxsamples_user_13 + Rxsamples_user_14 + Rxsamples_user_15 + Rxsamples_user_16 + Rxsamples_user_17 + Rxsamples_user_18 + Rxsamples_user_19 + Rxsamples_user_20 + Rxsamples_user_21 + Rxsamples_user_22 + Rxsamples_user_23 + Rxsamples_user_24 + noise;   %adding noise
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
            estimated_channel_matrix_u7 = special_impulse_search_function(M,N,l_p_u7,k_p_u7,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u8 = special_impulse_search_function(M,N,l_p_u8,k_p_u8,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u9 = special_impulse_search_function(M,N,l_p_u9,k_p_u9,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u10 = special_impulse_search_function(M,N,l_p_u10,k_p_u10,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u11 = special_impulse_search_function(M,N,l_p_u11,k_p_u11,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u12 = special_impulse_search_function(M,N,l_p_u12,k_p_u12,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u13 = special_impulse_search_function(M,N,l_p_u13,k_p_u13,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u14 = special_impulse_search_function(M,N,l_p_u14,k_p_u14,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u15 = special_impulse_search_function(M,N,l_p_u15,k_p_u15,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u16 = special_impulse_search_function(M,N,l_p_u16,k_p_u16,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u17 = special_impulse_search_function(M,N,l_p_u17,k_p_u17,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u18 = special_impulse_search_function(M,N,l_p_u18,k_p_u18,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u19 = special_impulse_search_function(M,N,l_p_u19,k_p_u19,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u20 = special_impulse_search_function(M,N,l_p_u20,k_p_u20,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u21 = special_impulse_search_function(M,N,l_p_u21,k_p_u21,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u22 = special_impulse_search_function(M,N,l_p_u22,k_p_u22,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u23 = special_impulse_search_function(M,N,l_p_u23,k_p_u23,l_max,k_max,threshold(th),SNR(k),Y_DD);
            estimated_channel_matrix_u24 = special_impulse_search_function(M,N,l_p_u24,k_p_u24,l_max,k_max,threshold(th),SNR(k),Y_DD);



            %MSE computation
            mse_impulse_pilot_user_1(th,k) = mse_impulse_pilot_user_1(th,k)+(((norm(estimated_channel_matrix_u1-true_channel_matrix_user_1,'fro'))^2)/(norm(true_channel_matrix_user_1,'fro')^2));
            mse_impulse_pilot_user_2(th,k) = mse_impulse_pilot_user_2(th,k)+(((norm(estimated_channel_matrix_u2-true_channel_matrix_user_2,'fro'))^2)/(norm(true_channel_matrix_user_2,'fro')^2));
            mse_impulse_pilot_user_3(th,k) = mse_impulse_pilot_user_3(th,k)+(((norm(estimated_channel_matrix_u3-true_channel_matrix_user_3,'fro'))^2)/(norm(true_channel_matrix_user_3,'fro')^2));
            mse_impulse_pilot_user_4(th,k) = mse_impulse_pilot_user_4(th,k)+(((norm(estimated_channel_matrix_u4-true_channel_matrix_user_4,'fro'))^2)/(norm(true_channel_matrix_user_4,'fro')^2));
            mse_impulse_pilot_user_5(th,k) = mse_impulse_pilot_user_5(th,k)+(((norm(estimated_channel_matrix_u5-true_channel_matrix_user_5,'fro'))^2)/(norm(true_channel_matrix_user_5,'fro')^2));
            mse_impulse_pilot_user_6(th,k) = mse_impulse_pilot_user_6(th,k)+(((norm(estimated_channel_matrix_u6-true_channel_matrix_user_6,'fro'))^2)/(norm(true_channel_matrix_user_6,'fro')^2));
            mse_impulse_pilot_user_7(th,k) = mse_impulse_pilot_user_7(th,k)+(((norm(estimated_channel_matrix_u7-true_channel_matrix_user_7,'fro'))^2)/(norm(true_channel_matrix_user_7,'fro')^2));
            mse_impulse_pilot_user_8(th,k) = mse_impulse_pilot_user_8(th,k)+(((norm(estimated_channel_matrix_u8-true_channel_matrix_user_8,'fro'))^2)/(norm(true_channel_matrix_user_8,'fro')^2));
            mse_impulse_pilot_user_9(th,k) = mse_impulse_pilot_user_9(th,k)+(((norm(estimated_channel_matrix_u9-true_channel_matrix_user_9,'fro'))^2)/(norm(true_channel_matrix_user_9,'fro')^2));        
            mse_impulse_pilot_user_10(th,k) = mse_impulse_pilot_user_10(th,k)+(((norm(estimated_channel_matrix_u10-true_channel_matrix_user_10,'fro'))^2)/(norm(true_channel_matrix_user_10,'fro')^2));
            mse_impulse_pilot_user_11(th,k) = mse_impulse_pilot_user_11(th,k)+(((norm(estimated_channel_matrix_u11-true_channel_matrix_user_11,'fro'))^2)/(norm(true_channel_matrix_user_11,'fro')^2));
            mse_impulse_pilot_user_12(th,k) = mse_impulse_pilot_user_12(th,k)+(((norm(estimated_channel_matrix_u12-true_channel_matrix_user_12,'fro'))^2)/(norm(true_channel_matrix_user_12,'fro')^2));
            mse_impulse_pilot_user_13(th,k) = mse_impulse_pilot_user_13(th,k)+(((norm(estimated_channel_matrix_u13-true_channel_matrix_user_13,'fro'))^2)/(norm(true_channel_matrix_user_13,'fro')^2));
            mse_impulse_pilot_user_14(th,k) = mse_impulse_pilot_user_14(th,k)+(((norm(estimated_channel_matrix_u14-true_channel_matrix_user_14,'fro'))^2)/(norm(true_channel_matrix_user_14,'fro')^2));
            mse_impulse_pilot_user_15(th,k) = mse_impulse_pilot_user_15(th,k)+(((norm(estimated_channel_matrix_u15-true_channel_matrix_user_15,'fro'))^2)/(norm(true_channel_matrix_user_15,'fro')^2));
            mse_impulse_pilot_user_16(th,k) = mse_impulse_pilot_user_16(th,k)+(((norm(estimated_channel_matrix_u16-true_channel_matrix_user_16,'fro'))^2)/(norm(true_channel_matrix_user_16,'fro')^2));
            mse_impulse_pilot_user_17(th,k) = mse_impulse_pilot_user_17(th,k)+(((norm(estimated_channel_matrix_u17-true_channel_matrix_user_17,'fro'))^2)/(norm(true_channel_matrix_user_17,'fro')^2));
            mse_impulse_pilot_user_18(th,k) = mse_impulse_pilot_user_18(th,k)+(((norm(estimated_channel_matrix_u18-true_channel_matrix_user_18,'fro'))^2)/(norm(true_channel_matrix_user_18,'fro')^2));
            mse_impulse_pilot_user_19(th,k) = mse_impulse_pilot_user_19(th,k)+(((norm(estimated_channel_matrix_u19-true_channel_matrix_user_19,'fro'))^2)/(norm(true_channel_matrix_user_19,'fro')^2));        
            mse_impulse_pilot_user_20(th,k) = mse_impulse_pilot_user_20(th,k)+(((norm(estimated_channel_matrix_u20-true_channel_matrix_user_20,'fro'))^2)/(norm(true_channel_matrix_user_20,'fro')^2));
            mse_impulse_pilot_user_21(th,k) = mse_impulse_pilot_user_21(th,k)+(((norm(estimated_channel_matrix_u21-true_channel_matrix_user_21,'fro'))^2)/(norm(true_channel_matrix_user_21,'fro')^2));
            mse_impulse_pilot_user_22(th,k) = mse_impulse_pilot_user_22(th,k)+(((norm(estimated_channel_matrix_u22-true_channel_matrix_user_22,'fro'))^2)/(norm(true_channel_matrix_user_22,'fro')^2));
            mse_impulse_pilot_user_23(th,k) = mse_impulse_pilot_user_23(th,k)+(((norm(estimated_channel_matrix_u23-true_channel_matrix_user_23,'fro'))^2)/(norm(true_channel_matrix_user_23,'fro')^2));        
            mse_impulse_pilot_user_24(th,k) = mse_impulse_pilot_user_24(th,k)+(((norm(estimated_channel_matrix_u24-true_channel_matrix_user_24,'fro'))^2)/(norm(true_channel_matrix_user_24,'fro')^2));

        end
    end
    end


mse_impulse_pilot_user_1 = mse_impulse_pilot_user_1/iter;
mse_impulse_pilot_user_2 = mse_impulse_pilot_user_2/iter;
mse_impulse_pilot_user_3 = mse_impulse_pilot_user_3/iter;
mse_impulse_pilot_user_4 = mse_impulse_pilot_user_4/iter;
mse_impulse_pilot_user_5 = mse_impulse_pilot_user_5/iter;
mse_impulse_pilot_user_6 = mse_impulse_pilot_user_6/iter;
mse_impulse_pilot_user_7 = mse_impulse_pilot_user_7/iter;
mse_impulse_pilot_user_8 = mse_impulse_pilot_user_8/iter;
mse_impulse_pilot_user_9 = mse_impulse_pilot_user_9/iter;
mse_impulse_pilot_user_10 = mse_impulse_pilot_user_10/iter;
mse_impulse_pilot_user_11 = mse_impulse_pilot_user_11/iter;
mse_impulse_pilot_user_12 = mse_impulse_pilot_user_12/iter;
mse_impulse_pilot_user_13 = mse_impulse_pilot_user_13/iter;
mse_impulse_pilot_user_14 = mse_impulse_pilot_user_14/iter;
mse_impulse_pilot_user_15 = mse_impulse_pilot_user_15/iter;
mse_impulse_pilot_user_16 = mse_impulse_pilot_user_16/iter;
mse_impulse_pilot_user_17 = mse_impulse_pilot_user_17/iter;
mse_impulse_pilot_user_18 = mse_impulse_pilot_user_18/iter;
mse_impulse_pilot_user_19 = mse_impulse_pilot_user_19/iter;
mse_impulse_pilot_user_20 = mse_impulse_pilot_user_20/iter;
mse_impulse_pilot_user_21 = mse_impulse_pilot_user_21/iter;
mse_impulse_pilot_user_22 = mse_impulse_pilot_user_22/iter;
mse_impulse_pilot_user_23 = mse_impulse_pilot_user_23/iter;
mse_impulse_pilot_user_24 = mse_impulse_pilot_user_24/iter;


figure;
semilogy(SNR_dB,mse_impulse_pilot_user_1(1,:),'r--o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_2(1,:),'r--+','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_3(1,:),'r--*','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_4(1,:),'r--.','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_5(1,:),'r--x','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_6(1,:),'r--_','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_7(1,:),'r--square','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_8(1,:),'r--diamond','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_9(1,:),'r--^','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_10(1,:),'r--v','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_11(1,:),'r--o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_12(1,:),'r--+','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_13(1,:),'r--*','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_14(1,:),'r--.','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_15(1,:),'r--x','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_16(1,:),'r--_','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_17(1,:),'r--square','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_18(1,:),'r--diamond','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_19(1,:),'r--^','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_20(1,:),'r--v','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_21(1,:),'r--square','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_22(1,:),'r--diamond','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_23(1,:),'r--^','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 
semilogy(SNR_dB,mse_impulse_pilot_user_24(1,:),'r--v','LineWidth',3,'MarkerFaceColor','r','MarkerSize',9.0);
hold on; grid on; 

axis tight;
title('Pilot-Impulse based CSI estimation with threshold = 3');
xlabel('SNR(dB)'); ylabel('NMSE');
legend('User 1','User 2','User 3','User 4','User 5','User 6','User 7','User 8','User 9','User 10''User 11','User 12','User 13','User 14','User 15','User 16','User 17','User 18','User 19','User 20','User 21','User 22','User 23','User 24');