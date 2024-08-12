function modified_with_loop_omp_estimation_vector = omp_with_side_channel_knowledge_estimation(received_vector,sensing_matrix,nonz)

M = size(sensing_matrix,1);
N = size(sensing_matrix,2);
modified_with_loop_omp_estimation_vector = zeros(N,1)+1j*zeros(N,1);

r_curr = received_vector;
% r_prev = zeros(M,1)+1j*zeros(M,1);
Q_matrix = [];
I_indices = [];
% iteration_count = 0;


while nonz > 0
%     iteration_count = iteration_count +1;
    [val,ind] = max(abs(sensing_matrix'*r_curr));   %correlation step
    I_indices = [I_indices ind];    %index finding for the most correlated column
    Q_matrix = [Q_matrix sensing_matrix(:,ind)];  %support matrix update
%     ls_est = (inv(Q_matrix'*Q_matrix))*(Q_matrix')*received_vector;
    ls_est = pinv(Q_matrix)*received_vector;   %perfect for time domain
%     placement
%     ls_est = Q_matrix\received_vector;
%     r_prev = r_curr;
    r_curr = received_vector - Q_matrix*ls_est;  %residue computation
    nonz = nonz-1;
end

modified_with_loop_omp_estimation_vector(I_indices) = ls_est;   %assigning to the final sparse vector




