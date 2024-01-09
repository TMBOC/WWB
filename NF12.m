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


%% plot Fig 12
K = K_vec(1);
SNR_vec = 10.^(SNR_dB_vec/10);

clc;close; figure(12);
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

x_text_vec = [-14, -9, -7, -18, -14];
y_text_vec = [3.5, 3, -1, -11, -13];
for idx_kappa=1:length(kappa_vec)
    hw_WWB=plot(SNR_dB_vec,WWB(idx_kappa,:,1,1,1),'Marker', '.','Linewidth',lw, 'color', C{idx_kappa},'MarkerSize',12);hold on;
    text(x_text_vec(idx_kappa), y_text_vec(idx_kappa), ['\kappa=' num2str(kappa_vec(idx_kappa))],'color', C{idx_kappa},'Fontsize',FS);hold on;
    hw_EZZB=plot(SNR_dB_vec,10*log10(sqrt(EZZB(idx_kappa,:))), '-','Linewidth',lw, 'color', C{idx_kappa});hold on;
end

line([-15, -14],[1.778,3],'color', C{1});
line([-14, -14],[0.4413, 3],'color', C{1});
line([-9, -11],[2,-0.8871],'color', C{2});
line([-9, -10],[2,-2.971],'color', C{2});
line([-7, -9],[-2,-3.801],'color', C{3});
line([-7, -8],[-2,-6.306],'color', C{3});
line([-16, -15],[-10.5,-4.162],'color', C{4});
line([-16, -14],[-10.5,-5.361],'color', C{4});
line([-11.7, -11],[-12.5,-9.098],'color', C{5});
line([-11.7, -10],[-12.5,-10.24],'color', C{5});


% % % % 
% % % % ann_x = [0.57 0.49];
% % % % ann_y = [0.55 0.515];
% % % % annotation('textarrow',ann_x,ann_y,'String','(-4, -10.25)','color', C{1},'Fontsize',FS);
% % % % ann_x = [0.57 0.495];
% % % % ann_y = [0.46 0.48];
% % % % annotation('line',ann_x,ann_y,'color', C{1});

xlabel('SNR (dB)','Fontsize',FS)
ylabel('10*log_{10}(RMSE)','Fontsize',FS)
legend([hw_WWB hw_EZZB, hw_BCRB],'WWB', 'ZZB','BCRB','Fontsize',FS,'Box','off')
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

save("WWB_NF11.mat")