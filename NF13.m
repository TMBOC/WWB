close all; clear; clc;


% von Mises
mu_vec = pi*0;
% kappa_vec = [0, 1, 2, 5, 20];
kappa_vec = [1];

% signal & noise sample
SNR_dB_vec = -20:1:15;
K_vec = 16;

group_vec = [10];


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% WWB
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


WWB = WWB_evaluate(kappa_vec, SNR_dB_vec, mu_vec, K_vec, group_vec);


%% ML & MAP
K=K_vec(1);k=0:(K-1);
theta_a = pi/6;

%******************
% Beta distribution
% a=20;
% von Mises
mu = mu_vec(1);
kappa = kappa_vec(1);
%******************

SNR_vec = 10.^(SNR_dB_vec/10);
sigma = 1;
 

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
EZZB = zeros(1, length(SNR_dB_vec));

L = 2^15; % grid size of fft
omega_d = 2*pi*((-L/2):(L/2-1))/L;
omega_da = circ_vmrnd(mu, kappa, n_trials)';

estError_ML=zeros(length(SNR_vec), n_trials);
estError_MAP=zeros(length(SNR_vec), n_trials);
for idx_SNR = 1:length(SNR_vec)
    h = waitbar( 0, sprintf( 'Monte-Carlo simulation for %.2f dB in progress ...', SNR_dB_vec(idx_SNR) ) );
    waitbar( idx_SNR/length(SNR_vec), h );
    
    omega_d_ML = zeros(1, length(n_trials));
    omega_d_MAP = zeros(1, length(n_trials));
%     omega_da = zeros(1, length(n_trials));
    
    b = sqrt(2 * SNR_vec(idx_SNR));

    
    parfor idb = 1:n_trials'
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
    close( h );
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



%% plot
close;clc;figure(13);
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

y_text_vec = [4, 2.6, 0.9, -2, -5.5];
hw_WWB=plot(SNR_dB_vec,WWB(1,:,1,1),... % kappa=1, all SNRs, K=16
                                    '-','Linewidth',lw, 'color', C{3},'marker', '.','MarkerSize',10);hold on;

BCRB_vm = zeros(1, length(SNR_vec));

for idx_SNR = 1:length(SNR_vec)
    SNR = SNR_vec(idx_SNR);
    J_F = SNR* (K*(K-1)*(2*K-1)) / 3;
%         ECRB = 1./J_F;
    BCRB_vm(idx_SNR) = 10*log10(sqrt(1./(J_F + kappa * besseli(1,kappa)/besseli(0,kappa))));
end
hw_BCRB = plot(SNR_dB_vec, BCRB_vm, '^-', 'color', C{5}); hold on;

% % % %     J_omega = (SNR)*K*(K-1)*(2*K-1)/3;
% % % %     sigma2 = pi*pi/(2*a(n)+1);
% % % %     ZZB_omega(n,:)=2*sigma2*(1-normcdf(sqrt(K*SNR),0,1))+(J_omega.^(-1)).*gammainc(K*SNR/2,1.5);


for idx_SNR = 1:length(SNR_vec)
    SNR = SNR_vec(idx_SNR);
    J_F = SNR* (K*(K-1)*(2*K-1)) / 3;
    sigma2 = -2* log(besseli(1,kappa)/besseli(0,kappa));
    EZZB(idx_SNR)=2*sigma2*(1-normcdf(sqrt(K*SNR),0,1))+(J_F.^(-1)).*gammainc(K*SNR/2,1.5);
end



%hw_ML = plot(SNR_dB_vec, 10*log10(RMSE_ML),'marker',"diamond", 'color', C{3}); hold on;
hw_MAP = plot(SNR_dB_vec, 10*log10(RMSE_MAP),'marker','*', 'color', C{1}); hold on;
hw_EZZB = plot(SNR_dB_vec, 10*log10(sqrt(EZZB)));
hold off

xlabel('SNR (dB)','Fontsize',FS)
ylabel('10*log_{10}(RMSE)','Fontsize',FS)
legendtext = {'WWB', 'ZZB', 'MAP', 'BCRB'};
legend([hw_WWB, hw_EZZB, hw_MAP, hw_BCRB],legendtext,'Fontsize',FS,'box','off');
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid on
axis([min(SNR_dB_vec) max(SNR_dB_vec) -25 5]);
ax = gca;
ax.LineWidth = 1.5;

save("WWB_NF13.mat")
%%
[A1, B1] = max(WWB-10*log10(sqrt(EZZB)))
