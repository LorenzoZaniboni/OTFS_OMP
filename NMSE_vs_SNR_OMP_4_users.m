close all;
clear all;
clc;
addpath(pwd + "/utils");



SNR_dB = -10:5:40;
SNR = 10.^(SNR_dB/10);
M = 32;
N = 32;
l_max = 8;   % max delay tap      %upto 12.5 microseconds for conventional system
k_max = 4;   % max doppler tap    %upto 500 KMPH doppler for conventional system
Ncp = l_max;
taps = 4;
iter = 100;
n_users = 4;
n_pilots = M*N;

mse_omp_time_placed_pilots = zeros(1,length(SNR_dB));
mse_omp_time_placed_pilots_with_side_channel_knowledge = zeros(1,length(SNR_dB));


%% 4 USER CASE
pilot_vector_1 = 2*randi([0,1],256,1)-1;      %bpsk pilot symbols for user 1
pilot_vector_2 = 2*randi([0,1],256,1)-1;      %bpsk pilot symbols for user 2
pilot_vector_3 = 2*randi([0,1],256,1)-1;      %bpsk pilot symbols for user 3
pilot_vector_4 = 2*randi([0,1],256,1)-1;      %bpsk pilot symbols for user 4

% Plot Tx Vector
pilot_vector_1 = [pilot_vector_1;zeros(768,1)];   
pilot_vector_2 = [zeros(256,1);pilot_vector_2;zeros(512,1)];   
pilot_vector_3 = [zeros(512,1);pilot_vector_3;zeros(256,1)];    
pilot_vector_4 = [zeros(768,1);pilot_vector_4];                 




%%
sensing_matrix_1 = sensing_matrix_creation(M,N,l_max,k_max,n_pilots,pilot_vector_1);
sensing_matrix_2 = sensing_matrix_creation(M,N,l_max,k_max,n_pilots,pilot_vector_2);
sensing_matrix_3 = sensing_matrix_creation(M,N,l_max,k_max,n_pilots,pilot_vector_3);
sensing_matrix_4 = sensing_matrix_creation(M,N,l_max,k_max,n_pilots,pilot_vector_4);
sensing_matrix = [sensing_matrix_1 sensing_matrix_2 sensing_matrix_3 sensing_matrix_4];

for ii = 1:iter

    delay_taps_1 = randperm(l_max,taps)-1;     %to get zero as a possible tap
    delay_taps_2 = randperm(l_max,taps)-1;
    delay_taps_3 = randperm(l_max,taps)-1;
    delay_taps_4 = randperm(l_max,taps)-1;
    doppler_taps_1 = randperm(k_max,taps)-1;   %to get zero as a possible tap   
    doppler_taps_2 = randperm(k_max,taps)-1;
    doppler_taps_3 = randperm(k_max,taps)-1;
    doppler_taps_4 = randperm(k_max,taps)-1;
    chann_coeff_1 = (1/sqrt(2))*(randn(1,taps) + 1j*randn(1,taps));
    chann_coeff_2 = (1/sqrt(2))*(randn(1,taps) + 1j*randn(1,taps));
    chann_coeff_3 = (1/sqrt(2))*(randn(1,taps) + 1j*randn(1,taps));
    chann_coeff_4 = (1/sqrt(2))*(randn(1,taps) + 1j*randn(1,taps));
    noise = (sqrt(1/2))*(randn(1,n_pilots) + 1j*randn(1,n_pilots));
    
    for k = 1:length(SNR_dB)
        disp("# of iteration: " + ii + "; SNR: " + SNR_dB(k))
        TxSamples_1 = sqrt(SNR(k))*pilot_vector_1.';
        TxSamples_2 = sqrt(SNR(k))*pilot_vector_2.';
        TxSamples_3 = sqrt(SNR(k))*pilot_vector_3.';
        TxSamples_4 = sqrt(SNR(k))*pilot_vector_4.';

        TxSamplesCP_1 = [TxSamples_1(n_pilots-Ncp+1:end) TxSamples_1];
        TxSamplesCP_2 = [TxSamples_2(n_pilots-Ncp+1:end) TxSamples_2];
        TxSamplesCP_3 = [TxSamples_3(n_pilots-Ncp+1:end) TxSamples_3];
        TxSamplesCP_4 = [TxSamples_4(n_pilots-Ncp+1:end) TxSamples_4];

        RxsamplesCP_1 = zeros(1,n_pilots+Ncp)+1j*zeros(1,n_pilots+Ncp);
        RxsamplesCP_2 = zeros(1,n_pilots+Ncp)+1j*zeros(1,n_pilots+Ncp);
        RxsamplesCP_3 = zeros(1,n_pilots+Ncp)+1j*zeros(1,n_pilots+Ncp);
        RxsamplesCP_4 = zeros(1,n_pilots+Ncp)+1j*zeros(1,n_pilots+Ncp);

        % Rx OTFS signals with CP
        for tap_iter = 1:taps
            RxsamplesCP_1 = RxsamplesCP_1 + chann_coeff_1(tap_iter)*circshift(TxSamplesCP_1.*exp(1j*2*pi/M*...
                (-Ncp:n_pilots-1)*doppler_taps_1(tap_iter)/N),[0,delay_taps_1(tap_iter)]);
            RxsamplesCP_2 = RxsamplesCP_2 + chann_coeff_2(tap_iter)*circshift(TxSamplesCP_2.*exp(1j*2*pi/M*...
                (-Ncp:n_pilots-1)*doppler_taps_2(tap_iter)/N),[0,delay_taps_2(tap_iter)]);
            RxsamplesCP_3 = RxsamplesCP_3 + chann_coeff_3(tap_iter)*circshift(TxSamplesCP_3.*exp(1j*2*pi/M*...
                (-Ncp:n_pilots-1)*doppler_taps_3(tap_iter)/N),[0,delay_taps_3(tap_iter)]);
            RxsamplesCP_4 = RxsamplesCP_4 + chann_coeff_4(tap_iter)*circshift(TxSamplesCP_4.*exp(1j*2*pi/M*...
                (-Ncp:n_pilots-1)*doppler_taps_4(tap_iter)/N),[0,delay_taps_4(tap_iter)]);

        end
        Rxsamples_1 = (RxsamplesCP_1(Ncp+1:n_pilots+Ncp));
        Rxsamples_2 = (RxsamplesCP_2(Ncp+1:n_pilots+Ncp));
        Rxsamples_3 = (RxsamplesCP_3(Ncp+1:n_pilots+Ncp));
        Rxsamples_4 = (RxsamplesCP_4(Ncp+1:n_pilots+Ncp));

        Rxsamples = (Rxsamples_1 + Rxsamples_2 + Rxsamples_3 + Rxsamples_4 + noise).'; 
        sensing_matrix_snr = sqrt(SNR(k)) * sensing_matrix;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Estimation Techniques
        
        h_sparse_omp = OMP_estimation(Rxsamples,sensing_matrix_snr,1);
        h_sparse_omp_with_side_channel_knowledge = omp_with_side_channel_knowledge_estimation(Rxsamples,sensing_matrix_snr,(n_users*taps));        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        h_sparse_actual_1 = h_sparse_creation(M,N,l_max,k_max,n_pilots,chann_coeff_1,taps,delay_taps_1,doppler_taps_1);
        h_sparse_actual_2 = h_sparse_creation(M,N,l_max,k_max,n_pilots,chann_coeff_2,taps,delay_taps_2,doppler_taps_2);
        h_sparse_actual_3 = h_sparse_creation(M,N,l_max,k_max,n_pilots,chann_coeff_3,taps,delay_taps_3,doppler_taps_3);
        h_sparse_actual_4 = h_sparse_creation(M,N,l_max,k_max,n_pilots,chann_coeff_4,taps,delay_taps_4,doppler_taps_4);

        h_sparse_actual = [h_sparse_actual_1;h_sparse_actual_2;h_sparse_actual_3;h_sparse_actual_4];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Analysis by MSE computation
        
        mse_omp_time_placed_pilots(k) = mse_omp_time_placed_pilots(k)+((norm(h_sparse_actual-h_sparse_omp))^2)/(norm(h_sparse_actual)^2);
        mse_omp_time_placed_pilots_with_side_channel_knowledge(k) = mse_omp_time_placed_pilots_with_side_channel_knowledge(k)+((norm(h_sparse_actual-h_sparse_omp_with_side_channel_knowledge))^2)/(norm(h_sparse_actual)^2);
    end
end

mse_omp_time_placed_pilots = mse_omp_time_placed_pilots/iter;
mse_omp_time_placed_pilots_with_side_channel_knowledge = mse_omp_time_placed_pilots_with_side_channel_knowledge/iter;


semilogy(SNR_dB,mse_omp_time_placed_pilots,'r-o','LineWidth',2,'MarkerFaceColor','r','MarkerSize',1.4);
hold on; 
semilogy(SNR_dB,mse_omp_time_placed_pilots_with_side_channel_knowledge,'k-s','LineWidth',2,'MarkerFaceColor','k','MarkerSize',1.4);
grid on;
axis tight;
legend('OMP','OMP with side channel knowledge');
xlabel('SNR [dB]');
ylabel('NMSE');
