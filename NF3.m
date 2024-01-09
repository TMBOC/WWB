% example of von Mises Distribution for eqn. (12)
close all; clc; clear;

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

figure(301);
mu_vec = 0;
theta_vec = [-1:0.01:1]*pi;
kappa_vec = [0, 1, 2, 5, 20];

for idx_kappa = 1:length(kappa_vec)
    mu = mu_vec(1);
    kappa = kappa_vec(idx_kappa);
    p_vec = 1/(2*pi*besseli(0,kappa)) * exp(kappa*cos(theta_vec-mu));
%     hw1(idx_kappa) = plot(theta_vec, p,'Linewidth',lw); hold on;
    hw1(idx_kappa) = plot(theta_vec, p_vec, 'Linewidth',lwW, 'color', C{idx_kappa}); hold on;
    legendtext{idx_kappa} = ['\kappa = ' num2str(kappa_vec(idx_kappa))];
end
hold off;
% xlabel('\theta / \pi','Fontsize',FS)
ylabel('Probability Density','Fontsize',FS)
legend(hw1,legendtext,'Fontsize',FS,'box','off');
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid on
axis([-pi pi 0.0 1.0]);
ax = gca;
ax.LineWidth = 1.5;
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ax.FontSize = FS;
%%
clc; clear;

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

figure(302);
mu_vec = [0, 1, 0.5, 0, -0.5]*pi;
theta_vec = [-1:0.01:1]*pi;
kappa_vec = [0, 1, 2, 5, 20];

legendnumberText = {'\kappa = 0,  \mu = 0', ...
                    '\kappa = 1,  \mu = \pi',...
                    '\kappa = 2,  \mu = \pi/2', ...
                    '\kappa = 5,  \mu = 0', ...
                    '\kappa = 20, \mu = - \pi/2'};
for idx_mu = 1:length(mu_vec)
    mu = mu_vec(idx_mu);
    kappa = kappa_vec(idx_mu);
    p_vec = 1/(2*pi*besseli(0,kappa)) * exp(kappa*cos(theta_vec-mu));
%     hw1(idx_kappa) = plot(theta_vec, p,'Linewidth',lw); hold on;
    hw1(idx_mu) = plot(theta_vec, p_vec, 'Linewidth',lwW, 'color', C{idx_mu}); hold on;
%     legendtext{idx_mu} = legendnumberText;
end
hold off;
% xlabel('\theta / \pi','Fontsize',FS)
ylabel('Probability Density','Fontsize',FS)
% legend(hw1,legendnumberText,'Fontsize',FS,'box','off');
legend(hw1,legendnumberText,'Fontsize',10,'box','off');
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid on
axis([-pi pi 0.0 1.0]);
ax = gca;
ax.LineWidth = 1.5;
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ax.FontSize = FS;
