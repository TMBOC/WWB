close all; clear; clc;


% von Mises
mu_vec = pi*(-1:0.2:1);
kappa_vec = 0:1:10;

% signal & noise sample
SNR_dB_vec = [0 -10];
SNR_vec = 10.^(SNR_dB_vec/10); % SNR
K_vec = [100, 200];

% WWB setting
s_vec = 0.5;
s = s_vec(1);
si=s;
sj=s;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% WWB
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dd = 0.0001; % for integration
WWB = zeros(length(K_vec), length(SNR_vec), ...
            length(s_vec));% 3rd dimension for 's'. Note we set si=sj here.


for idx_K_SNR = 1:length(K_vec)
    K = K_vec(idx_K_SNR);
    SNR = SNR_vec(idx_K_SNR);

    h = waitbar( 0, sprintf( 'WWB evaluation for K = %d, SNR = %d dB,\n %s = %d %s rad and %s = %d rad^-2', K_vec(1), SNR_dB_vec(1),'\mu', round(mu_vec(1)/pi), '\pi','\kappa',kappa_vec(1) ) );
    cont=1/length(mu_vec)/length(kappa_vec);
    for idx_mu = 1:length(mu_vec)
        mu = mu_vec(idx_mu);
        for idx_kappa = 1:length(kappa_vec)
            waitbar( cont, h, ...
                sprintf( 'WWB evaluation for K = %d, SNR = %d dB,\n %s = %d %s rad and %s = %d rad^{-2}', K_vec(idx_K_SNR), SNR_dB_vec(idx_K_SNR), '\mu',round(mu_vec(idx_mu)/pi), '\pi','\kappa',kappa_vec(idx_kappa) ));
%                 sprintf( 'WWB evaluation for K = %d, SNR = %d dB, mu = %f rad and kappa = %d rad^2', K_vec(idx_K), SNR_dBvec(idx_SNR), mu_vec(idx_mu), kappa_vec(idx_kappa) ));
            cont = cont+1/length(mu_vec)/length(kappa_vec);

            kappa = kappa_vec(idx_kappa);

            % Setting up Test points
            SL = zeros(1,K/2-1);% test points ('S'), Sidelobe peaks
            for k=1:K/2-1
                SL(k) = 2*(k+0.5-0.25*(1-k/(K/2-1)))/K;
            end
            rp = 0.1:0.1:1; % test points ('E'), Evenly-distributed
            hv = [0.001 0.01 SL rp]*pi; % test points 'C'+'S'+'E'

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
            
%             Q = (Q+Q')/2;
            WWB(idx_mu, idx_kappa, idx_K_SNR) = 10*log10(sqrt(hv*inv(Q)*hv.'));
        end % kappa
    end % mu
    close(h);
end % (K,SNR) pair


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

[X,Y] = meshgrid(mu_vec, kappa_vec);
figure(901);
surf(X',Y', WWB(:,:,1));
xlim([-pi pi]);
shading(gca, 'interp');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\mu (rad)', 'FontSize', FS)
ylabel('\kappa (rad^{-2})', 'FontSize', FS)
zlabel('10*log_{10}(RMSE)', 'FontSize', FS)
ax = gca;
ax.FontSize = 14; 
az = 135;
el = 45;
view(az, el);

figure(902);
surf(X',Y', WWB(:,:,2));
shading(gca, 'interp');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
xlabel('\mu (rad)', 'FontSize', FS)
xlim([-pi pi]);
ylabel('\kappa (rad^{-2})', 'FontSize', FS)
zlabel('10*log_{10}(RMSE)', 'FontSize', FS)
ax = gca;
ax.FontSize = 14; 
az = 135;
el = 45;
view(az, el);

save("WWB_NF9.mat")