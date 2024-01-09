% WWB adapted from NF5.m
close all; clear; clc;


% von Mises
mu = pi/2;
kappa = 2;

% signal & noise sample
SNR_dB_vec = -20 : 1 : 10;
SNR_vec = 10.^(SNR_dB_vec/10); % SNR
K_vec = [20];
K=K_vec(1);

% WWB settings
s_vec = [0.5];
s = s_vec(1);
si=s;
sj=s;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% WWB
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dd = 0.0001; % for integration


% Setting up Test points
SL = zeros(1,K/2-1);% test points ('S'), Sidelobe peaks
for k=1:K/2-1
    SL(k) = 2*(k+0.5-0.25*(1-k/(K/2-1)))/K;
end
hh = [0.001 0.01 SL]*pi; % test points 'C'+'S'; Our purpose is just to see differences between number of test points; using 'E' elements will not introduce much insight
nh = length(hh);
WWB = zeros(nh-2, length(SNR_vec));% excluding 0.001*pi, 0.01*pi, temporarily, we will add them back during evaluation


h = waitbar( 0, sprintf( 'WWB evaluation for  SNR = %d dB and number of test points = 2+%d ', SNR_dB_vec(1),1) );
cont=1/(nh-2)/length(SNR_vec)*2;
for idx_nh = 1:2:(nh-2)
    fprintf('testing %d ''S'' elements\n', idx_nh);
    for idx_SNR = 1:length(SNR_vec)
        waitbar(cont, h, sprintf( 'WWB evaluation for  SNR = %d dB and number of test points = 2 + %d ', SNR_dB_vec(idx_SNR),idx_nh));
        cont = cont + 1/(nh-2)/length(SNR_vec)*2;
        SNR = SNR_vec(idx_SNR);

        hv = hh(1:(idx_nh+2)); % include 0.001*pi, 0.01*pi and [1, 3, 5, 7, 9]

        Q = zeros(length(hv),length(hv));
        for ii=1:length(hv)
            hi = hv(ii);

            % denominator, part 1 of 2
            mu_hi = -si*(1-si)*2*SNR*(K-cos(hi*(K-1)/2)*sin(hi*K/2)/sin(hi/2));
            d = hi/(2*pi);
            d = round(d/dd)*dd;
            v = (0+dd/2): dd : (1-d-dd/2);
            ph = exp(kappa * ((si-1)*cos(2*pi*v-mu) - si*cos(2*pi*v-mu+hi))) / besseli(0, kappa);
            exp_gam_hi = sum(ph)*dd;
            exp_eta_hi = exp(mu_hi) * exp_gam_hi;
            for jj=1:length(hv)
                hj = hv(jj);

                % denominator, part 2 of 2
                mu_hj = -sj*(1-sj)*2*SNR*(K-cos(hj*(K-1)/2)*sin(hj*K/2)/sin(hj/2));
                d = hj/(2*pi);
                d = round(d/dd)*dd;
                v = (0+dd/2): dd : (1-d-dd/2);
                ph = exp(kappa * ((sj-1)*cos(2*pi*v-mu) - sj*cos(2*pi*v-mu+hj))) / besseli(0, kappa);
                exp_gam_hj = sum(ph)*dd;
                exp_eta_hj = exp(mu_hj) * exp_gam_hj;

                % nominator
                if jj==ii
                    mu_4 = SNR*(...
                                K*((si+sj-1)^2+(si-1)^2+(sj-1)^2-1) ...
                                -2*(si+sj-1)*(si-1)*cos(hi*(K-1)/2)*sin(hi*K/2)/sin(hi/2) ...
                                -2*(si+sj-1)*(sj-1)*cos(hj*(K-1)/2)*sin(hj*K/2)/sin(hj/2) ...
                                +2*(si-1)*(sj-1)*K...
                                );
                    mu_1 = SNR*(...
                                K*((si+sj-1)^2+si^2+sj^2-1) ...
                                +2*si*sj*K...
                                -2*(si+sj-1)*si*cos(hi*(K-1)/2)*sin(hi*K/2)/sin(hi/2) ...
                                -2*(si+sj-1)*sj*cos(hj*(K-1)/2)*sin(hj*K/2)/sin(hj/2)...
                                );
                else    
                    mu_4 = SNR*(...
                                K*((si+sj-1)^2+(si-1)^2+(sj-1)^2-1) ...
                                -2*(si+sj-1)*(si-1)*cos(hi*(K-1)/2)*sin(hi*K/2)/sin(hi/2) ...
                                -2*(si+sj-1)*(sj-1)*cos(hj*(K-1)/2)*sin(hj*K/2)/sin(hj/2) ...
                                +2*(si-1)*(sj-1)*cos((hi-hj)*(K-1)/2)*sin((hi-hj)*K/2)/sin((hi-hj)/2)...
                                );
                    mu_1 = SNR*(...
                                K*((si+sj-1)^2+si^2+sj^2-1) ...
                                +2*si*sj*cos((hi-hj)*(K-1)/2)*sin((hi-hj)*K/2)/sin((hi-hj)/2)...
                                -2*(si+sj-1)*si*cos(hi*(K-1)/2)*sin(hi*K/2)/sin(hi/2) ...
                                -2*(si+sj-1)*sj*cos(hj*(K-1)/2)*sin(hj*K/2)/sin(hj/2)...
                                );
                end
                mu_2 = SNR*(...
                            K*(sj^2+(si-1)^2+(si-sj)^2-1) ...
                            -2*sj*(si-1)*cos((hi+hj)*(K-1)/2)*sin((hi+hj)*K/2)/sin((hi+hj)/2)...
                            +2*sj*(si-sj)*cos(hj*(K-1)/2)*sin(hj*K/2)/sin(hj/2)...
                            -2*(si-1)*(si-sj)*cos(hi*(K-1)/2)*sin(hi*K/2)/sin(hi/2)...
                            );
                mu_3 = SNR*(...
                            K*(si^2+(sj-1)^2+(si-sj)^2-1) ...
                            -2*si*(sj-1)*cos((hi+hj)*(K-1)/2)*sin((hi+hj)*K/2)/sin((hi+hj)/2)...
                            +2*si*(sj-si)*cos(hi*(K-1)/2)*sin(hi*K/2)/sin(hi/2)...
                            -2*(sj-1)*(sj-si)*cos(hj*(K-1)/2)*sin(hj*K/2)/sin(hj/2)...
                            );

                d = abs(hi-hj)/2/pi; 
                d = round(d/dd)*dd;
                d4 = min(hi,hj)/2/pi;
                v = (d4+dd/2) : dd : (1-d-dd/2);
                ph4 = (1/besseli(0, kappa)) * exp(kappa * ...
                                                        ((1-si-sj)*cos(2*pi*(v+d)-mu) ...
                                                          + (si-1)*cos(2*pi*(v+d)-mu-hi) ...
                                                          + (sj-1)*cos(2*pi*(v+d)-mu-hj)));
                exp_gam_4 = sum(ph4)*dd;

                d = abs(hi-hj)/2/pi;
                d = round(d/dd)*dd;
                d1 = min(hi,hj)/2/pi;
                v = (d1+dd/2) : dd : (1-d-dd/2);
                ph1 = (1/besseli(0, kappa)) * exp(kappa * ...
                                                        ((si+sj-1)*cos(2*pi*v-mu-hj) ...
                                                          - si*cos(2*pi*(v+d)-mu) ...
                                                          - sj*cos(2*pi*v-mu))...
                        );
                exp_gam_1 = sum(ph1)*dd;

                d = (hi+hj)/2/pi;
                d = round(d/dd)*dd;
                d2 = 0;
                v = (d2+dd/2) : dd : (1-d-dd/2);
                ph2 = (1/besseli(0, kappa)) * exp(kappa * ...
                                                        ((sj-si)*cos(2*pi*v+hi-mu) ...
                                                          - sj*cos(2*pi*(v+d)-mu)...
                                                          + (si-1)*cos(2*pi*v-mu)));
                exp_gam_2 = sum(ph2)*dd;

                d = (hi+hj)/2/pi;
                d = round(d/dd)*dd;
                d3 = 0;
                v = (d3+dd/2) : dd : (1-d-dd/2);
                ph3 = (1/besseli(0, kappa)) * exp(kappa * ...
                                                        ((si-sj)*cos(2*pi*v+hj-mu) ...
                                                          - si*cos(2*pi*(v+d)-mu) ...
                                                          + (sj-1)*cos(2*pi*v-mu)));
                exp_gam_3 = sum(ph3)*dd;

                exp_eta_1 = exp(mu_1) * exp_gam_1;
                exp_eta_2 = exp(mu_2) * exp_gam_2;
                exp_eta_3 = exp(mu_3) * exp_gam_3;
                exp_eta_4 = exp(mu_4) * exp_gam_4;

                Q(ii,jj) = (exp_eta_1 - exp_eta_2 - exp_eta_3 + exp_eta_4)/(exp_eta_hi*exp_eta_hj);
            end % jj
        end % ii
        WWB(idx_nh, idx_SNR) = 10*log10(sqrt(hv*inv(Q)*hv.'));
    end % SNR
end % number of test points, i.e. length of subsets of  hv
close(h);

%% plot Fig 8
close; clc;figure(806); % figure 8(f)
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

Lt = {'-', '-.','--', '*-', '^-'}; % Marker and line types

idx_plot = 1;
for idx_nh=1:2:nh-2
    hw1(idx_plot)=plot(SNR_dB_vec,WWB(idx_nh,:),Lt{idx_plot},'Linewidth',lw, 'color', C{idx_plot});hold on;
    legendText{idx_plot} = [num2str(idx_nh) ' points'];
    idx_plot = idx_plot + 1;
end
hold off
xlabel('SNR (dB)','Fontsize',FS)
ylabel('10*log_{10}(RMSE)','Fontsize',FS)
ax=gca;
ax.FontSize = 14; 

legend(hw1, legendText,'Box','off', 'Fontsize', FS);
set(gca,'XMinorTick','on','YMinorTick','on');
grid on
axis([min(SNR_dB_vec) max(SNR_dB_vec) -25 5]);
save("WWB_NF7f.mat")