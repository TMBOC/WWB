function [WWB] = WWB_evaluate(kappa_vec, SNR_dB_vec, mu_vec, K_vec, group_vec)
%[Functions] evaluate Weiss-Weinstein Bound (WWB) for digitial frequency
%estimation considering von Mises distribution as prior distribution.
%
%[Dependencies] this function requires only MATLAB built-in functions
%
%WWB = WWB_evaluate(SNR_vec, kappa_vec, mu_vec, group_vec)
%
%   Inputs:    
%       SNR_vec        - 1xS vector of SNR in dB
%       kappa_vec      - 1xV vector of the concentration parameter of the
%                        von Mises distribution
%       mu_vec         - 1xV vector of the mean parameter of the von Mises 
%                        distribution. Default value is 0.
%       K_vec          - 1xK vector of how many samples in one period of
%                        the baseband signal. Default value is 16.
%       group_vec      - 1xG vector of numbers of 'E' elements in P evaluations of the WWB.
%                        example: [10 0] means that two groups of WWB will
%                        be evaluated: one of them is (2,9,10) while the
%                        other is (2,9,0). For the moment we only support
%                        2 groups: [numOfE, 0]
%   Outputs:
%       WWB            - a 3d matrix of dimension VxSxG storing WWB for (idx_kappa, idx_SNR, idx_group)
%                        where idx_kappa, idx_SNR, idx_group are index for
%                        kappa_vec, SNR_vec, and group_vec
% Author: Xin Zhang 2023-12-27 
%--------------------------------------------------------------------------
% WWB setting
s_vec = 0.5;
s = s_vec(1);
si=s;
sj=s;
dd = 0.00001; % for integration

WWB = zeros(length(kappa_vec), length(SNR_dB_vec), length(K_vec), length(group_vec));% 3rd dimension for comparing (C,S,E) and (C,S,0).   

SNR_vec = 10.^(SNR_dB_vec/10); % SNR

h = waitbar( 0, sprintf( 'WWB evaluation for group %d\n SNR = %d dB and %s = %d rad^{-2}',1, SNR_dB_vec(1),'\kappa',kappa_vec(1) ) );
cont=1/length(kappa_vec)/length(SNR_vec)/length(group_vec)/length(K_vec);

for idx_group = 1:1:length(group_vec) 
    if group_vec(idx_group)~=0
        flag_group = 1; % group 1 (C,S,E)
    else
        flag_group = 2; % group 2(C,S,0)
    end
    for idx_K = 1:length(K_vec)
        K = K_vec(idx_K);
        for idx_SNR = 1:length(SNR_vec)
            SNR = SNR_vec(idx_SNR);
            for idx_kappa = 1:length(kappa_vec)% kappa and mu are always set simultaneously
                kappa = kappa_vec(idx_kappa);
                mu = mu_vec(idx_kappa); % kappa and mu are always set simultaneously

                waitbar(cont, h, sprintf( 'WWB evaluation for group %d\n SNR = %d dB, %s = %d rad^{-2}, %s = %1.1f, K = %d', idx_group, SNR_dB_vec(idx_SNR),'\kappa',kappa_vec(idx_kappa),'\mu',mu_vec(idx_kappa), K_vec(idx_K) ));
                cont = cont+1/length(kappa_vec)/length(SNR_vec)/length(group_vec)/length(K_vec);

                % Setting up Test points
                SL = zeros(1,K/2-1);% test points ('S'), Sidelobe peaks
                for k=1:K/2-1
                    SL(k) = 2*(k+0.5-0.25*(1-k/(K/2-1)))/K;
                end
                
                if (flag_group == 1)
                    resolution_rp = round(1/group_vec(idx_group), 1);
                    rp = 0.1:resolution_rp:1; % test points ('E'), Evenly-distributed
                    fprintf('group_vec(%d) = %d. We have %d elements ''E'' with resolution_rp = %1.1f within [0.1, 1]%s \n',idx_group,group_vec(idx_group),length(rp),resolution_rp,'\pi');
                    hv = [0.001 0.01 SL rp]*pi; % test points 'C'+'S'+'E'
                    
                else
                    fprintf('group_vec(%d) = %d. No elements ''E'' \n',idx_group,group_vec(idx_group));
                    hv = [0.001 0.01 SL]*pi; % test points 'C'+'S'
                end

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

%                 Q = (Q+Q')/2;
                WWB(idx_kappa, idx_SNR, idx_K, idx_group) = 10*log10(sqrt(hv*inv(Q)*hv.'));
            end % kappa
        end % SNR
    end % K
end % group
close(h);