% Adapted from NF8.m
close all; clear; clc;


% von Mises

kappa_vec = [0, 1, 2, 5, 20];
mu_vec = pi*zeros(1,length(kappa_vec));

% signal & noise sample
SNR_dB_vec = -20:1:15;
K_vec = 20;

group_vec = [10 0];


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% WWB
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


WWB = WWB_evaluate(kappa_vec, SNR_dB_vec, mu_vec, K_vec, group_vec);

%% plot Fig 7(a)
close; figure(701);
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

y_text_vec = [4, 2.6, 0.9, -2, -5.5]-1;
idx_K = 1;
for idx_kappa=1:length(kappa_vec)
    hw1=plot(SNR_dB_vec,WWB(idx_kappa,:,idx_K,1),'Marker', '.','Linewidth',lw-1, 'color', C{idx_kappa},'MarkerSize',10);
    text(-19, y_text_vec(idx_kappa), ['\kappa=' num2str(kappa_vec(idx_kappa))],'color', C{idx_kappa},'Fontsize',FS);hold on;
    hw2=plot(SNR_dB_vec,WWB(idx_kappa,:,idx_K,2), '-','Linewidth',lw, 'color', C{idx_kappa});
end
hold off

ann_x = [0.57 0.49];
ann_y = [0.55 0.515];
annotation('textarrow',ann_x,ann_y,'String','(-4, -10.25)','color', C{1},'Fontsize',FS);
ann_x = [0.57 0.495];
ann_y = [0.46 0.48];
annotation('textarrow',ann_x,ann_y,'String','(-4, -11.2)','color', C{1},'Fontsize',FS);

xlabel('SNR (dB)','Fontsize',FS)
ylabel('10*log_{10}(RMSE)','Fontsize',FS)
legend([hw1 hw2],'Proposed: WWB_{(2,9,10)}', 'Legacy: WWB_{(2,9,0)}','Fontsize',FS-2,'box','off')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out','Fontsize',FS)
grid on
axis([min(SNR_dB_vec) max(SNR_dB_vec) -25 5]);
ax = gca;
ax.LineWidth = 1.5;

%% plot Fig 7(b)
close;figure(702);
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
x_text_vec_zoom = [-4.9, -6, -6.9, -5.95, -7];
y_text_vec_zoom = [-7.2, -7.2, -7.2, -10.4, -12.2];
idx_K = 1;
for idx_kappa=1:length(kappa_vec)
    hw1=plot(SNR_dB_vec,WWB(idx_kappa,:,idx_K,1),'Marker', '.','Linewidth',lw, 'color', C{idx_kappa},'MarkerSize',15);
    text(x_text_vec_zoom(idx_kappa), y_text_vec_zoom(idx_kappa), ['\kappa=' num2str(kappa_vec(idx_kappa))],'color', C{idx_kappa},'Fontsize',FS);hold on;
    hw2=plot(SNR_dB_vec,WWB(idx_kappa,:,idx_K,2), '-','Linewidth',lw, 'color', C{idx_kappa});
end
hold off

ann_x = [0.57 0.525]+0.055;
ann_y = [0.855 0.855]-0.31;
annotation('textarrow',ann_x,ann_y,'String','(-4, -10.25)','color', C{1},'Fontsize',FS);
ann_x = [0.57 0.525]+0.05;
ann_y = [0.765 0.765]-0.33;
annotation('textarrow',ann_x,ann_y,'String','(-4, -11.2)','color', C{1},'Fontsize',FS);

xlabel('SNR (dB)','Fontsize',FS)
ylabel('10*log_{10}(RMSE)','Fontsize',FS)
legend([hw1 hw2],'Proposed: WWB_{(2,9,10)}', 'Legacy: WWB_{(2,9,0)}','Fontsize',FS, 'box', 'off')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out','Fontsize',FS)
grid on
ax = gca;
ax.LineWidth = 1.5;

axis([-8 -1 -14 -7]);

%% 
% Maximum difference, index of SNR_dB_vec
[A1, B1] = max(abs(WWB(:,:,idx_K,1)-WWB(:,:,idx_K,2)));
[A2, B2] = max(A1) % print it!

save("WWB_NF7.mat")
