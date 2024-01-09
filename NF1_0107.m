%% Figure 1
clear; clc;
K_vec=16;
K=K_vec(1);
k=0:(K-1);
theta_a = pi/6;

%******************
% Beta distribution
% a=20;
% von Mises
mu_vec = 0;
mu = mu_vec(1);
kappa_vec = 20;
kappa = kappa_vec(1);
%******************

SNR_dB_vec = -20:15;
SNR_vec = 10.^(SNR_dB_vec/10);
sigma = 1;
group_vec=[10];
WWB = WWB_evaluate(kappa_vec, SNR_dB_vec, mu_vec, K_vec, group_vec);
%%
J_D = SNR_vec * (K*(K-1)*(2*K-1)) / 3;
ECRB = 1./J_D;
%******************
% Beta distribution
% BCRB_beta = 1./(J_D + (a-1)*(2*a-1)/(pi^2*(a-2)));
% von Mises
BCRB_vm = 1./(J_D + kappa * besseli(1,kappa)/besseli(0,kappa));
%******************

n_trials = 10000;
RMSE_ML = zeros(1, length(SNR_dB_vec));
RMSE_MAP = zeros(1, length(SNR_dB_vec));

L = 2^15; % grid size of fft
omega_d = 2*pi*((-L/2):(L/2-1))/L;
omega_da = circ_vmrnd(mu, kappa, n_trials)';

estError_ML=zeros(length(SNR_vec), n_trials);
estError_MAP=zeros(length(SNR_vec), n_trials);
h = waitbar( 0, sprintf( 'Monte-Carlo simulation for %.2f dB in progress ...', SNR_dB_vec(1) ) );
cont = 1/length(SNR_vec);
for idx_SNR = 1:length(SNR_vec)
    waitbar( cont, h, sprintf( 'Monte-Carlo simulation for %.2f dB in progress ...', SNR_dB_vec(idx_SNR) ) );
    cont = cont + 1/length(SNR_vec);
    
    omega_d_ML = zeros(1, length(n_trials));
    omega_d_MAP = zeros(1, length(n_trials));
%     omega_da = zeros(1, length(n_trials));
    
    b = sqrt(2 * SNR_vec(idx_SNR));

    
    parfor idb = 1:n_trials
% % % %         % my implementation of generating omega_da
% % % %         t = random('Beta', a, a);
% % % %         omega_da(idb) = 2*pi*t - pi;
% % % % %         % reference implementation of generating omega_da
% % % % %         y = rand;
% % % % %         omega_da(idb) = 2*pi*betainv(y, a, a) - pi;
        s_tilde = b * exp(1i*(omega_da(idb)*k + theta_a));
        w_tilde = sigma * (randn(1, K) + 1i * randn(1,K));
        x_tilde = s_tilde + w_tilde;
        F = fft(x_tilde, L)/K;
        F = fftshift(F);
        AS = 2*K*SNR_vec(idx_SNR)*real(exp(-1i * theta_a) * F);
        [~, idmax_ML] = max(AS);
        omega_d_ML(idb) = omega_d(idmax_ML);
        A_2 = AS + ...
                kappa*cos(omega_d - mu);
%               (a-1)*log((pi+omega_d)/(2*pi)) + (a-1)*log((pi-omega_d)/(2*pi));
        [~, idmax_MAP] = max(A_2);
        omega_d_MAP(idb) = omega_d(idmax_MAP);
    end
    
    estError_ML(idx_SNR, :) = omega_d_ML - omega_da;
    estError_MAP(idx_SNR, :) = omega_d_MAP - omega_da;
    for idx_estError1 = 1:n_trials
        if (estError_ML(idx_SNR, idx_estError1) > pi)
            estError_ML(idx_SNR, idx_estError1) = estError_ML(idx_SNR, idx_estError1) - 2*pi;
        end
        if (estError_ML(idx_SNR, idx_estError1) < -pi)
            estError_ML(idx_SNR, idx_estError1) = estError_ML(idx_SNR, idx_estError1) + 2*pi;
        end        
    end
    for idx_estError2 = 1:n_trials
        if (estError_MAP(idx_SNR, idx_estError2) > pi)
            estError_MAP(idx_SNR, idx_estError2) = estError_MAP(idx_SNR, idx_estError2) - 2*pi;
        end
        if (estError_MAP(idx_SNR, idx_estError2) < -pi)
            estError_MAP(idx_SNR, idx_estError2) = estError_MAP(idx_SNR, idx_estError2) + 2*pi;
        end        
    end
        
    RMSE_ML(idx_SNR) = sqrt((1/n_trials) * sum(estError_ML(idx_SNR, :).^2));
    RMSE_MAP(idx_SNR) = sqrt((1/n_trials) * sum(estError_MAP(idx_SNR, :).^2));
%     RMSE_ML(idx_SNR) = sqrt((1/n_trials) * (2-2*sum(cos(omega_d_ML - omega_da))));
%     RMSE_MAP(idx_SNR) = sqrt((1/n_trials) * (2-2*sum(cos(omega_d_MAP - omega_da))));
end
close( h );
save("WWB_NF1.mat")
%%
presentation = 1;      % 0 for manuscript (B/W), 1 for pres (color)
if presentation
    lw = 1.5;
    lwW = 2.0;
    FS = 14;
else
    lw = 0.5;
    lwW = 0.5;
    FS = 12;
end

C = {[0.3010 0.7450 0.9330],...
    'k',...
    [0.4660 0.6740 0.1880],...
    [0.4940 0.1840 0.5560],...
    [0.9290 0.6940 0.1250]}; % Cell array of colors.

close;figure(1);

plot(SNR_dB_vec, 10*log10(RMSE_ML), '-.', 'color', [0.6350 0.0780 0.1840], 'Linewidth', lwW);hold on;
plot(SNR_dB_vec, 10*log10(RMSE_MAP), '--', 'color', [0.6350 0.0780 0.1840], 'Linewidth', lwW);hold on;
plot(SNR_dB_vec, 10*log10(sqrt(ECRB)), 'marker','.','MarkerSize',14,'Color', C{3},'Linewidth', lwW);hold on;
plot(SNR_dB_vec, 10*log10(sqrt(BCRB_vm)), '-','Color',C{5},'Linewidth', lwW);hold on;
plot(SNR_dB_vec, WWB, '-', 'color', C{1},  'Linewidth', lwW);hold off;
grid on;axis([-20, 15, -25, 5]);
legend('ML','MAP','CRB','BCRB','WWB','Fontsize',10,'box','off');% CRB here is actually ECRB
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
xlabel('SNR (dB)');ylabel('10*log_{10}(RMSE)');
ax = gca;
ax.LineWidth = 1.5;
% ann_x = [0.25 0.15];
% ann_y = [0.4 0.57];
% annotation('line',ann_x,ann_y,'color', 'k');
% dim = [.1 .56 .07 .07];
% annotation('ellipse',dim)

% % Zoom-in
% axes('Position',[.2 .2 .2 .2])
% box on
% hw_zoom_ECRB=plot(SNR_dB_vec, 10*log10(sqrt(ECRB)), 'marker','.','MarkerSize',15,'Color', C{3});hold on;
% hw_zoom_BCRB=plot(SNR_dB_vec, 10*log10(sqrt(BCRB_vm)), '-','Color',C{5});
% grid on;axis([-20, -18, -7.5, -6.5]);
% legend('ECRB','BCRB');
% hold off;

