function [h_est,delay_taps_est,Doppler_taps_est,num_est_taps,estimated_channel_matrix] = Pilot_impulse_CE_modified(M,N,k_p,l_p,l_max,k_max,thrld,Y_DD,pp_dB)
    estimated_channel_matrix = zeros(l_max,k_max);
    h_est = [];
    delay_taps_est = [];
    Doppler_taps_est = [];
    num_est_taps = 0;
    for r = k_p:k_p+l_max-1
        for c = l_p:l_p+k_max-1
             r_mod = mod(r,M); %%%
             c_mod = mod(c,N);  %%%
              if r_mod == 0
                  r_mod = M;
              end
              if c_mod == 0
                  c_mod = N;
              end
              
              if(abs(Y_DD(r_mod,c_mod)) > thrld)
%                 estimated_channel_matrix(r-l_p+1,c-k_p+1) = Y_DD(r_mod,c_mod)/sqrt(10^(pp_dB/10));
                if r_mod - 1 < (r-k_p)
                    estimated_channel_matrix(r-k_p+1,c-l_p+1) = Y_DD(r_mod,c_mod)/((sqrt(10^(pp_dB/10)))*exp(-1i*2*pi*(c_mod-1)/N)*exp(1i*2*pi*(c-l_p)*mod((r_mod - 1 - (r-k_p)),M)/(M*N)));
                else 
                    estimated_channel_matrix(r-k_p+1,c-l_p+1) = Y_DD(r_mod,c_mod)/((sqrt(10^(pp_dB/10)))*exp(1i*2*pi*(c-l_p)*mod((r_mod -1 - (r-k_p)),M)/(M*N)));
                end
                %h_est = [h_est Y_DD(r,c)/sqrt(10^(pp_dB/10))];
                h_est = [h_est estimated_channel_matrix(r-k_p+1,c-l_p+1)];
                delay_taps_est = [delay_taps_est r-k_p];
                Doppler_taps_est = [Doppler_taps_est c-l_p];
                num_est_taps = num_est_taps + 1;
             end
        end
    end
end