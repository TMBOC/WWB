%% Figure 5， adapted from Figure 4(a)
close all; clc; clear;
figure(5)% normalized signal ambiguity surface, known phase, K=16, omega_da=0.4pi, theta_a=pi/6


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



b = 1; % normalized
omega_da = 0 * pi;
theta_a = pi/6;
K = 20;
% omega_d = pi * (-1:0.001:1);
L = 2^18;                   % FFT grid size
omega_d = 2*pi*[-(L/2):1:(L/2)-1]/L;
d_omega = omega_d-omega_da;
A_sn = (b/K) * cos(d_omega*(K-1)/2) .* sin(d_omega*K/2) ./ sin(d_omega/2);
plot(omega_d/pi, A_sn, 'Linewidth',lw);hold on;
xlabel('\theta/\pi','FontSize',FS); ylabel('A_{X}(\theta)','FontSize',FS);
axis([-1, 1, -0.2, 1]);grid on;



% Test points(part 1): Center point
C_x = 0.001*5; %5 is multiplied as a tolerance to better visualize a point with abscissa other than zero
C_y = (b/K) * cos((C_x*pi - omega_da)*(K-1)/2) .* sin((C_x*pi - omega_da)*K/2) ./ sin((C_x*pi - omega_da)/2);
C_x_1 = 0.01; %0.01 is added as a tolerance to better visualize a point with abscissa other than zero
C_y_1 = (b/K) * cos((C_x_1*pi - omega_da)*(K-1)/2) .* sin((C_x_1*pi - omega_da)*K/2) ./ sin((C_x_1*pi - omega_da)/2);
% plot(C_x,C_y,'or');
plot(C_x,C_y,'-s','MarkerSize',5,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);hold on;
plot(C_x_1,C_y_1,'-s','MarkerSize',5,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);hold on;
% x = [0.6 0.53];
% y = [0.88 0.88];
% a=annotation('textarrow',x,y,'String','''C''');
% a.Color='red';

% Test points(part 2): Sidelobe points
for k=1:K/2-1
    S_x(k) = 2*(k+0.5-0.25*(1-k/(K/2-1)))/K;
end
d_omega_new = S_x * pi - omega_da;
S_y = (b/K) * cos(d_omega_new*(K-1)/2) .* sin(d_omega_new*K/2) ./ sin(d_omega_new/2);
plot(S_x,S_y,'v','MarkerSize',5,...
    'MarkerEdgeColor','magenta',...
    'MarkerFaceColor',[0.4940 0.1840 0.5560]);hold on;

% Test points(part 3): Evenly distributed points
E_x = 0.1:0.1:1;
E_y = zeros(1, length(E_x));
plot(E_x,E_y,'^','MarkerSize',5,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .8 .1]);hold on;

plot([0.6],[0.9],'-s','MarkerSize',5,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);hold on;
plot([0.6],[0.82],'v','MarkerSize',5,...
    'MarkerEdgeColor','magenta',...
    'MarkerFaceColor',[0.4940 0.1840 0.5560]);hold on;
plot([0.6],[0.74],'^','MarkerSize',5,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .8 .1]);hold on;
xt = [0.7 0.7  0.7 ];
yt = [0.9 0.82 0.74];
str = {'''C''','''S''','''E'''};
text(xt,yt,str, 'Fontsize',FS);


hold off;