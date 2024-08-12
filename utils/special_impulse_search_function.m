function estimated_channel_matrix = special_impulse_search_function(M,N,l_p,k_p,l_max,k_max,threshold,SNR,Y_DD)

estimated_channel_matrix = zeros(l_max,k_max);

for r_u = l_p:l_p+l_max-1
    for c_u = k_p:k_p+k_max-1
        r_mod = mod(r_u,M);
        c_mod = mod(c_u,N);
        if r_mod == 0
           r_mod = M;
        end
        if c_mod == 0
            c_mod = N;
        end
        if abs(Y_DD(r_mod,c_mod)) > threshold
            if r_mod - 1 < (r_u-l_p)
                estimated_channel_matrix(r_u-l_p+1,c_u-k_p+1) = Y_DD(r_mod,c_mod)/((sqrt(SNR))*exp(-1i*2*pi*(c_mod-1)/N)*exp(1i*2*pi*(c_u-k_p)*mod((r_mod - 1 - (r_u-l_p)),M)/(M*N)));
            else
                estimated_channel_matrix(r_u-l_p+1,c_u-k_p+1) = Y_DD(r_mod,c_mod)/((sqrt(SNR))*exp(1i*2*pi*(c_u-k_p)*mod((r_mod -1 - (r_u-l_p)),M)/(M*N)));
            end
        end
    end
end
end