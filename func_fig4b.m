function [A_x, idmax] = func_fig4b(SNR_dB, omega_d, omega_da, theta_a, K)

M=3; % number of Monte Carlo trials
SNR = 10^(SNR_dB/10);
k = 0:K-1;
L=length(omega_d);
b=sqrt(2*SNR);


idmax = zeros(M,1);
A_x = zeros(M, L);
for m = 1:M
    nn = randn(1,K)+1i*randn(1,K);
    x = b*exp(1i*(omega_da*k+theta_a))+nn;
    
    F = fft(x,L)/K;
    F = fftshift(F);
    AS = real(exp(-1i*theta_a)*F);     % known theta ambiguity surface
    [~,idmax(m)]=max(AS);
    A_x(m,:)=AS;
end


end