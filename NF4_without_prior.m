%% Figure 4(a)
close all; clc; clear;
figure(401)% normalized signal ambiguity surface, known phase, K=16, omega_da=0.4pi, theta_a=pi/6
b = 1; % normalized
omega_da = 0.4 * pi;
theta_a = pi/6;
K = 20;
% omega_d = pi * (-1:0.001:1);
L = 2^15;                   % FFT grid size
omega_d = 2*pi*[-(L/2):1:(L/2)-1]/L;
d_omega = omega_d-omega_da;
A_sn = (b/K) * cos(d_omega*(K-1)/2) .* sin(d_omega*K/2) ./ sin(d_omega/2);
plot(omega_d/pi, A_sn);
xlabel('\theta/\pi'); ylabel('A_{X}(\theta)');
axis([-1, 1, -0.2, 1]);grid on;


%% Figure 4(b)



% preparing for plots
clc;
for i=1:3
    legendtext{i} = ['Trial ' num2str(i)];
end
% parameters to be estimated
%omega_d = pi * (-1:0.001:1);%duplicated in the above cell


figure(402); % Noisy ambiguity surface, A_x(omega_d), known phase. K=16, omega_da=0.4pi, theta_a=pi/6
subplot(3,1,1);
SNR_dB = 5; % this parameter changes with subplots
[A_x, idmax] = func_fig4b(SNR_dB, omega_d, omega_da, theta_a, K);
pl1 = plot(omega_d/pi, A_x(1,:), '-.', ...
           omega_d/pi, A_x(2,:), '--', ...
           omega_d/pi, A_x(3,:), '-'); hold on;
plot(omega_d(idmax(1))/pi,A_x(1,idmax(1)), '*', ...
     omega_d(idmax(2))/pi,A_x(2,idmax(2)), '*', ...
     omega_d(idmax(3))/pi,A_x(3,idmax(3)), '*');hold off;
legend(pl1, legendtext); 
title('SNR = 5 dB');ylabel('A_x(\theta)');

subplot(3,1,2)
SNR_dB = -10; % this parameter changes with subplots
[A_x, idmax] = func_fig4b(SNR_dB, omega_d, omega_da, theta_a, K);
pl1 = plot(omega_d/pi, A_x(1,:), '-.', ...
           omega_d/pi, A_x(2,:), '--', ...
           omega_d/pi, A_x(3,:), '-'); hold on;
plot(omega_d(idmax(1))/pi,A_x(1,idmax(1)), '*', ...
     omega_d(idmax(2))/pi,A_x(2,idmax(2)), '*', ...
     omega_d(idmax(3))/pi,A_x(3,idmax(3)), '*');hold off;
title('SNR = -10 dB');ylabel('A_x(\theta)');

subplot(3,1,3)
SNR_dB = -20; % this parameter changes with subplots
[A_x, idmax] = func_fig4b(SNR_dB, omega_d, omega_da, theta_a, K);
pl1 = plot(omega_d/pi, A_x(1,:), '-.', ...
           omega_d/pi, A_x(2,:), '--', ...
           omega_d/pi, A_x(3,:), '-'); hold on;
plot(omega_d(idmax(1))/pi,A_x(1,idmax(1)), '*', ...
     omega_d(idmax(2))/pi,A_x(2,idmax(2)), '*', ...
     omega_d(idmax(3))/pi,A_x(3,idmax(3)), '*');hold off; 
title('SNR = -20 dB');xlabel('\theta/\pi');ylabel('A_x(\theta)');


