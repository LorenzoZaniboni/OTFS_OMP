function h_sparse_actual = h_sparse_creation(M,N,l_max,k_max,n_pilots,h,taps,delay_taps,doppler_taps)

h_sparse_actual = zeros(l_max*k_max,1);
index_vector = zeros(taps,1);

for tap = 1:taps
    index_vector(tap) = k_max*delay_taps(tap) + doppler_taps(tap) +1;
end

h_sparse_actual(index_vector) = h;







