close all; clear; clc
%% Figures to compare DD-LQR, Input-Constrained DD-LQR, and RH IC DD-LQR 
% Algorithm 2 in ref [2], Algorithm 2, and Algorithm 3, respectively

Q = 1;
R = 1;

%% load data
load ddlqr
X_lqr  = X_log;
U_lqr  = U_log;
TG_lqr = Gamma;

load ddlqr_ic
X_ic  = X_log;
U_ic  = U_log;
T_ic  = Tau_log;
TG_ic = Tau_Gamma_log;

load ddlqr_rh
X_rh  = X_log;
U_rh  = U_log;
T_rh  = Tau_log;
TG_rh = Tau_Gamma_log;

clear X_log U_log Tau_log Tau_Gamma_log Gamma

t = length(U_ic);

color1 = [0         0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];
color3 = [0.9290    0.6940    0.1250];

%% Figure 1: Trajectories
b = 6;

figure()
plot(X_lqr(1,:), X_lqr(2,:), '*--', 'LineWidth', 1, 'Color', color1)
hold on
% plot(X_ic(1,:),  X_ic(2,:),  '*--', 'LineWidth', 1, 'Color', color3)
plot(X_rh(1,:),  X_rh(2,:),  '*--', 'LineWidth', 1, 'Color', color2)
% legend('ddlqr','ddlqr-ic','ddlqr-rh', 'Location', 'southwest')
legend('ddlqr','ours-updating $\mathcal{C}$', 'Location', 'southeast', 'Interpreter', 'latex')
title('Trajectories')
xlabel('x_1')
ylabel('x_2')
grid on
axis square
xlim([-b, b])
ylim([-b, b])
hold off

ax = gca; % Get the current axes handle

ax.Title.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Legend.FontSize = 16;

%% Figure 2: Control 

figure()
plot(1:t, U_lqr,'*--', 'LineWidth', 1, 'Color', color1)
hold on
% plot(1:t, U_ic, '*--', 'LineWidth', 1, 'Color', color3)
plot(1:t, U_rh, '*--', 'LineWidth', 1, 'Color', color2)
% legend('ddlqr','ddlqr-ic','ddlqr-rh', 'Location', 'southeast')
legend('ddlqr','ours-updating $\mathcal{C}$', 'Location', 'southeast', 'Interpreter', 'latex')
title('Control')
xlabel('t')
ylabel('u')
grid on
hold off

ax = gca; % Get the current axes handle

ax.Title.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Legend.FontSize = 16;

%% Figure 3: Worst-case upper bound TG decrease

figure()
plot(1:t, TG_ic, '*--', 'LineWidth', 1, 'Color', color3)
hold on
plot(1:t, TG_rh, '*--', 'LineWidth', 1, 'Color', color2)
plot(1:t, TG_lqr * ones(1,t), '*--', 'LineWidth', 1, 'Color', color1)
legend('ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
    'ddlqr', 'Location', 'northeast', 'Interpreter', 'latex')
title('Worst-case upper bound decrease')
xlabel('t')
ylabel('upper bound')
grid on

ax = gca; % Get the current axes handle

ax.Title.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Legend.FontSize = 16;


%% Figure 4: Consistency set shrinkage
compare_consistency