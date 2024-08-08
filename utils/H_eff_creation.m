function H_eff = H_eff_creation(M,N,taps,chan_coef,delay_taps,Doppler_taps,Prx,Ptx)
    H_mat = zeros(M*N,M*N);     %time domain channel 
    F_N = 1/sqrt(N)*dftmtx(N);
    omega = exp(1j*2*pi/(M*N));
    for tap_iter = 1:taps
        H_mat = H_mat + chan_coef(tap_iter)*circshift(eye(M*N),delay_taps(tap_iter))*...
            (diag(omega.^((0:M*N-1)*Doppler_taps(tap_iter))));
    end
    H_eff = kron(F_N,Prx)*H_mat*kron(F_N',Ptx);   %delay doppler channel
end