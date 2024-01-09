%% Figure 3(a)
close all; clc; clear;
% figure(301)% normalized signal ambiguity surface, known phase, K=16, omega_da=0.4pi, theta_a=pi/6
% b = 1; % normalized
omega_da = 0.4 * pi;
omega_d = [-1:0.001:1]*pi;
theta_a = pi/6;
K = 16;
% omega_d = pi * (-1:0.001:1);
L = 2^15;                   % FFT grid size
% omega_d = 2*pi*[-(L/2):1:(L/2)-1]/L;
% d_omega = omega_d-omega_da;
% A_sn = (b/K) * cos(d_omega*(K-1)/2) .* sin(d_omega*K/2) ./ sin(d_omega/2);
% plot(omega_d/pi, A_sn);
% xlabel('\theta/\pi'); ylabel('A_{X}(\theta)');
% axis([-1, 1, -0.2, 1]);grid on;


presentation = 0;      % 0 for manuscript (B/W), 1 for pres (color)
if presentation
    lw = 1.5;
    FS = 14;
else
    lw = 0.5;
    FS = 12;
end

lc = ['b';'r';'g'];

kappa = 5;
% mu = -0.5*pi;
mu = 0*pi;

SNR_vec = 10.^([5 -10 -20]/10); % SNR
ns = length(SNR_vec);
LP = zeros(ns,length(omega_d));
figure(302);
for n=1:ns
    subplot(3,1,n)
    AS = 2*K*SNR_vec(n)*real(sum(exp(1i*[0:1:K-1].'*(omega_d-omega_da))))/K;
%     PS = kappa * (omega_d-mu) - log(2*pi*besseli(0,kappa));
    sigma2 = -2*log(besseli(1,kappa)/besseli(0,kappa));
    PS = -0.5*((omega_d-mu).^2./sigma2+log(2*pi*sigma2));
    LP(n,:) = AS+PS;
    a1 = ceil(max(LP(n,:)));
    a2 = floor(min(LP(n,:)));
    plot(omega_d/pi,LP(n,:),lc(n,:),'Linewidth',lw)
    grid
    axis([-1 1 -10 max(a1,10)])
    title(['SNR = ' num2str(10*log10(SNR_vec(n))) ' dB'],'Fontsize',FS)
    ylabel('A_{XP}(\theta)','Fontsize',FS)
    set(gca,'XTick',[-1:0.2:1])
    set(gca,'Fontsize',FS)
end
xlabel('\theta/\pi','Fontsize',FS)
