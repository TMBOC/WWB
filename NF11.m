% WWB and plots adapted from NF6.m, EZZB adapted from NF11.m
close all; clear; clc;


% von Mises

kappa_vec = [0, 1, 2, 5, 20];
mu_vec = pi*zeros(1,length(kappa_vec));

% signal & noise sample
SNR_dB_vec = -20:1:15;
K_vec = 20;


group_vec = 10;


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% WWB
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


WWB = WWB_evaluate(kappa_vec, SNR_dB_vec, mu_vec, K_vec, group_vec);



%% ML & MAP
clc;close all;
mv_vec = 0*pi;
mu = mu_vec(1);
K = K_vec(1);
k=0:(K-1);
theta_a = pi/6;
sigma = 1;
SNR_vec = 10.^(SNR_dB_vec/10);
n_trials = 10000;
RMSE_ML = zeros(length(kappa_vec), length(SNR_dB_vec));
RMSE_MAP = zeros(length(kappa_vec), length(SNR_dB_vec));

L = 2^15; % grid size of fft
omega_d = 2*pi*((-L/2):(L/2-1))/L;

h = waitbar( 0, sprintf( 'Monte-Carlo simulation for SNR = %d dB and %s = %d in progress ...', SNR_dB_vec(1) ,'\kappa', kappa_vec(1)));

cont = 1/length(kappa_vec)/length(SNR_vec);
for idx_kappa = 1:length(kappa_vec)
    kappa = kappa_vec(idx_kappa);
    omega_da = circ_vmrnd(mu, kappa, n_trials)';

    estError_ML=zeros(length(SNR_vec), n_trials);
    estError_MAP=zeros(length(SNR_vec), n_trials);
    for idx_SNR = 1:length(SNR_vec)
        
        waitbar( cont, h, sprintf( 'Monte-Carlo simulation for %.2f dB and %s = %d in progress ...', SNR_dB_vec(idx_SNR),'\kappa', kappa_vec(idx_kappa)));
        cont = cont + 1/length(kappa_vec)/length(SNR_vec);
        
        omega_d_ML = zeros(1, length(n_trials));
        omega_d_MAP = zeros(1, length(n_trials));
    %     omega_da = zeros(1, length(n_trials));

        b = sqrt(2 * SNR_vec(idx_SNR));


        parfor idb = 1:n_trials
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

        RMSE_ML(idx_kappa, idx_SNR) = sqrt((1/n_trials) * sum(estError_ML(idx_SNR, :).^2));
        RMSE_MAP(idx_kappa, idx_SNR) = sqrt((1/n_trials) * sum(estError_MAP(idx_SNR, :).^2));
    end % SNR_vec

end % kappa_vec

close(h);

%% plot Fig 10
K = K_vec(1);
SNR_vec = 10.^(SNR_dB_vec/10);

clc;close; figure(10);
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

EZZB = zeros(length(kappa_vec), length(SNR_dB_vec));
for idx_kappa=1:length(kappa_vec)
    for idx_SNR = 1:length(SNR_vec)
        SNR = SNR_vec(idx_SNR);
        J_F = SNR* (K*(K-1)*(2*K-1)) / 3;
%         sigma2 = 1 -
%         besseli(1,kappa_vec(idx_kappa))/besseli(0,kappa_vec(idx_kappa));% circular variance, improper here!!
        if kappa_vec(idx_kappa)==0
            sigma2 = 4*pi*pi/12;
        else
            sigma2 = -2 * log(besseli(1,kappa_vec(idx_kappa))/besseli(0,kappa_vec(idx_kappa))); % linear variance mentioned in circStat
        end
        EZZB(idx_kappa,idx_SNR)=2*sigma2*(1-normcdf(sqrt(K*SNR),0,1))+(J_F.^(-1)).*gammainc(K*SNR/2,1.5);
    end
end

BCRB_vm = zeros(length(kappa_vec), length(SNR_vec));
for idx_kappa = 1:length(kappa_vec)
    kappa = kappa_vec(idx_kappa);
    for idx_SNR = 1:length(SNR_vec)
        SNR = SNR_vec(idx_SNR);
        J_F = SNR* (K*(K-1)*(2*K-1)) / 3;
        BCRB_vm(idx_kappa, idx_SNR) = 10*log10(sqrt(1./(J_F + kappa * besseli(1,kappa)/besseli(0,kappa))));
    end
    hw_BCRB = plot(SNR_dB_vec, BCRB_vm(idx_kappa,:), '--', 'color', C{idx_kappa}); hold on;
end

x_text_vec = [-19, -19, -19, -19, -19]-0.5;
y_text_vec = [3.5, 1.7, 0.3, -2.6, -5.4];
for idx_kappa=1:length(kappa_vec)
    hw_WWB=plot(SNR_dB_vec,WWB(idx_kappa,:,1,1,1),'Marker', '.','Linewidth',lw, 'color', C{idx_kappa},'MarkerSize',12);hold on;
    text(x_text_vec(idx_kappa), y_text_vec(idx_kappa), ['\kappa=' num2str(kappa_vec(idx_kappa))],'color', C{idx_kappa},'Fontsize',12);hold on;
    hw_MAP = plot(SNR_dB_vec, 10*log10(RMSE_MAP(idx_kappa,:)),'marker','*', 'color', C{idx_kappa}); hold on;
end

% line([-15, -14],[2.242,3.5],'color', C{1});
% line([-13, -14],[1.297,3.5],'color', C{1});
% line([-10, -11],[2,0.7792],'color', C{2});
% line([-10, -8],[2,-2.2976],'color', C{2});
% line([-7, -9],[-2,-3.801],'color', C{3});
% line([-7, -8],[-2,-6.306],'color', C{3});
% line([-16, -15],[-10.5,-4.162],'color', C{4});
% line([-16, -14],[-10.5,-5.361],'color', C{4});
% line([-11.7, -11],[-12.5,-9.098],'color', C{5});
% line([-11.7, -10],[-12.5,-10.24],'color', C{5});


xlabel('SNR (dB)','Fontsize',FS)
ylabel('10*log_{10}(RMSE)','Fontsize',FS)
legend([hw_WWB,  hw_BCRB, hw_MAP],'WWB', 'BCRB','MAP','Fontsize',FS,'Box','off')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid on
axis([min(SNR_dB_vec) max(SNR_dB_vec) -25 5]);
ax = gca;
ax.LineWidth = 1.5;

ann_x = [0.25 0.14];
ann_y = [0.4 0.53];
annotation('line',ann_x,ann_y,'color', 'k');
dim = [.1 .52 .05 .05];
annotation('ellipse',dim)




% Zoom-in
axes('Position',[.2 .2 .2 .2])
box on
x_text_zoom_vec = [-19.5, -19.2, -18.9, -19.6, -19.9];
y_text_zoom_vec = [-8.4, -8.6, -8.8, -9.15, -9.4];
for idx_kappa = 1:length(kappa_vec)
    hw_zoom_BCRB=plot(SNR_dB_vec, BCRB_vm(idx_kappa, :), '--','Color',C{idx_kappa});hold on;
    text(x_text_zoom_vec(idx_kappa), y_text_zoom_vec(idx_kappa), ['\kappa=' num2str(kappa_vec(idx_kappa))],'color', C{idx_kappa});hold on;
end
grid on;axis([-20, -18, -9.5, -8]);
legend('BCRB', 'Box', 'off');

hold off

save("WWB_NF10.mat")