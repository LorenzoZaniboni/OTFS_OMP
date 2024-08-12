function sensing_matrix = sensing_matrix_creation(M,N,l_max,k_max,n_pilots,pilot_vector)

sensing_matrix = zeros(n_pilots,l_max*k_max);
delta_matrix = zeros(n_pilots,1);
identity_matrix= eye(n_pilots);
storage_matrix = zeros(n_pilots,n_pilots,l_max*k_max);

for l_i = 0:l_max-1
    if l_i == 0
        delta_matrix = exp((1j*2*pi/(M*N))*(0:n_pilots-1));
%         delta_matrix = diag(delta_matrix);
    else
        delta_matrix(1:n_pilots-l_i) = exp((1j*2*pi/(M*N))*(0:n_pilots-l_i-1));
        delta_matrix(n_pilots-l_i+1:end) = exp((1j*2*pi/(M*N))*(-l_i:-1));    %an n_pilots*1 vector
%         delta_matrix = diag(delta_matrix);   %we have delta_l_i matrix... an n_pilots*n_pilots matrix
    end
    for k_i = 0:k_max-1
%         delta_matrix = diag(delta_matrix.^k_i);    %raised to the power k_i
        cir_shifted_delta_matrix = circshift(identity_matrix,l_i)*diag(delta_matrix.^k_i);  %an n_pilots*n_pilots matrix
        storage_matrix(:,:,l_i*k_max+k_i+1) = cir_shifted_delta_matrix;
    end
end

for count = 1:l_max*k_max
    sensing_matrix(:,count) = storage_matrix(:,:,count)*pilot_vector;   %its a vector in the sening matrix
end







































end