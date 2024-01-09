function [A_x, idmax] = func_fig4b_vm(SNR_dB, omega_d, omega_da, theta_a, K, kappa, mu)

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
    AS = 2*K*SNR*real(exp(-1i*theta_a)*F);     % known theta ambiguity surface
%     PS = kappa * (omega_d-mu) - log(2*pi*besseli(0,kappa));
    PS = -0.5*((omega_d-mu).^2./sigma2+log(2*pi*sigma2))
    A_2 = AS+PS;
    [~,idmax(m)]=max(A_2);
    A_x(m,:)=AS;
end


end