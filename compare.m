close all; clear; %clc
%% Figures to compare DeePC, DD-LQR, Input-Constrained DD-LQR, and RH IC DD-LQR 
% DeePC, Algorithm 2 in ref [2], Algorithm 2, and Algorithm 3, respectively

load alldata Q R UU

%% load data
load mdlqr_ic
if size(UU,1) > 1
    U_log = vecnorm(UU*U_log, 2);
end
X_md = X_log;
U_md = U_log;
G_md = Gamma_log;
V_md = V_log;

t = length(U_md);

load deepc
if size(UU,1) > 1
    U_log = vecnorm(UU*U_log, 2);
end
X_dpc = X_log;
U_dpc = U_log;

load ddlqr
if size(UU,1) > 1
    U_log = vecnorm(UU*U_log, 2);
end
X_lqr = X_log;
U_lqr = U_log;
G_lqr = Gamma * ones(1, t);
V_lqr = V_log;

load ddlqr_ic
if size(UU,1) > 1
    U_log = vecnorm(UU*U_log, 2);
end
X_ic  = X_log;
U_ic  = U_log;
T_ic  = Tau_log;
G_ic  = Gamma_log;
V_ic  = V_log;

load ddlqr_rh
if size(UU,1) > 1
    U_log = vecnorm(UU*U_log, 2);
end
X_rh  = X_log;
U_rh  = U_log;
T_rh  = Tau_log;
G_rh  = Gamma_log;
V_rh  = V_log;

clear X_log U_log Tau_log Gamma_log Gamma V_log


color1 = [0         0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];
color3 = [0.9290    0.6940    0.1250];
color4 = [0.4940    0.1840    0.5560];
color5 = [0.4660    0.6740    0.1880];

%% Figure 1: Trajectories
b = 6;

figure()
plot(X_lqr(1,:), X_lqr(2,:), '*--', 'LineWidth', 1, 'Color', color1)
hold on
plot(X_md(1,:),  X_md(2,:),  '*--', 'LineWidth', 1, 'Color', color3)
plot(X_rh(1,:),  X_rh(2,:),  '*--', 'LineWidth', 1, 'Color', color2)
plot(X_dpc(1,:), X_dpc(2,:), '*--', 'LineWidth', 1, 'Color', color4)
% legend('ddlqr','ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
%     'deepc','Location', 'southeast', 'Interpreter', 'latex')
legend('ddlqr','ours model-based','ours data-driven',...
    'deepc','Location', 'southeast', 'Interpreter', 'latex')
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


%% Figure 1-1: Trajectories coordinate 34

if size(X_rh, 1) > 2

    b = 3;
    
    figure()
    plot(X_lqr(3,:), X_lqr(4,:), '*--', 'LineWidth', 1, 'Color', color1)
    hold on
    plot(X_md(3,:),  X_md(4,:),  '*--', 'LineWidth', 1, 'Color', color3)
    plot(X_rh(3,:),  X_rh(4,:),  '*--', 'LineWidth', 1, 'Color', color2)
    plot(X_dpc(3,:), X_dpc(4,:), '*--', 'LineWidth', 1, 'Color', color4)
    % legend('ddlqr','ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
    %     'deepc','Location', 'southeast', 'Interpreter', 'latex')
    legend('ddlqr','ours model-based','ours data-driven',...
        'deepc','Location', 'southeast', 'Interpreter', 'latex')
    title('Trajectories')
    xlabel('x_3')
    ylabel('x_4')
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

end

%% Figure 1-2: Trajectories
%{
b = 6;

T=t+1;
figure()
hold on
h11 =   plot(1:T, X_lqr(1,:), '*--', 'LineWidth', 1, 'Color', color1);
h12 =   plot(1:T, X_lqr(2,:), '*-',  'LineWidth', 1, 'Color', color1);
h13 =   plot(1:T, X_lqr(3,:), '*:',  'LineWidth', 1, 'Color', color1);
h14 =   plot(1:T, X_lqr(4,:), '*-.', 'LineWidth', 1, 'Color', color1);

h21 =   plot(1:T, X_ic(1,:),  '*--', 'LineWidth', 1, 'Color', color3);
        plot(1:T, X_ic(2,:),  '*-',  'LineWidth', 1, 'Color', color3);
        plot(1:T, X_ic(3,:),  '*:',  'LineWidth', 1, 'Color', color3);
        plot(1:T, X_ic(4,:),  '*-.', 'LineWidth', 1, 'Color', color3);

h31 =   plot(1:T, X_rh(1,:),  '*--', 'LineWidth', 1, 'Color', color2);
        plot(1:T, X_rh(2,:),  '*-',  'LineWidth', 1, 'Color', color2);
        plot(1:T, X_rh(3,:),  '*:',  'LineWidth', 1, 'Color', color2);
        plot(1:T, X_rh(4,:),  '*-.', 'LineWidth', 1, 'Color', color2);

h41 =   plot(1:T, X_dpc(1,:), '*--', 'LineWidth', 1, 'Color', color4);
        plot(1:T, X_dpc(2,:), '*-',  'LineWidth', 1, 'Color', color4);
        plot(1:T, X_dpc(3,:), '*:',  'LineWidth', 1, 'Color', color4);
        plot(1:T, X_dpc(4,:), '*-.', 'LineWidth', 1, 'Color', color4);


legend([h11,h21,h31,h41], {'ddlqr','ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$','deepc'},...
    'Location', 'northeast', 'Interpreter', 'latex')
title('Trajectories')
xlabel('t')
ylabel('x')
grid on
axis square
% xlim([-b, b])
% ylim([-b, b])
hold off

ax = gca; % Get the current axes handle

ax.Title.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Legend.FontSize = 16;

%}

%% Figure 2: Control 

figure()
plot(1:t, U_lqr,'*--', 'LineWidth', 1, 'Color', color1)
hold on
plot(1:t, U_md, '*--', 'LineWidth', 1, 'Color', color3)
plot(1:t, U_rh, '*--', 'LineWidth', 1, 'Color', color2)
plot(1:t, U_dpc,'*--', 'LineWidth', 1, 'Color', color4)
% legend('ddlqr','ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
%     'deepc', 'Location', 'northeast', 'Interpreter', 'latex')
legend('ddlqr','ours model-based','ours data-driven',...
    'deepc','Location', 'northeast', 'Interpreter', 'latex')
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
hold on
plot(1:t, G_lqr,'*--', 'LineWidth', 1, 'Color', color1)
plot(1:t, G_md, '*--', 'LineWidth', 1, 'Color', color3)
plot(1:t, G_rh, '*--', 'LineWidth', 1, 'Color', color2)
plot(1:t, G_ic, '*--', 'LineWidth', 1, 'Color', color5)
% legend('ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
%     'ddlqr', 'Location', 'east', 'Interpreter', 'latex')
legend('ddlqr','ours model-based','ours data-driven',...
    'ours nonRH data-driven','Location', 'east', 'Interpreter', 'latex')
title('Worst-case upper bound decrease')
xlabel('t')
ylabel('upper bound')
grid on

ax = gca; % Get the current axes handle

ax.Title.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Legend.FontSize = 16;


%% Figure 3-1: Lyapunov function V=x'*inv(P)*x decrease

figure()
hold on
plot(1:t, V_lqr,'*--', 'LineWidth', 1, 'Color', color1)
plot(1:t, V_md, '*--', 'LineWidth', 1, 'Color', color3)
plot(1:t, V_rh, '*--', 'LineWidth', 1, 'Color', color2)
% plot(1:t, V_ic, '*--', 'LineWidth', 1, 'Color', color5)
% legend('ddlqr','ours model-based','ours data-driven',...
%     'ours nonRH data-driven','Location', 'northeast', 'Interpreter', 'latex')
legend('ddlqr','ours model-based','ours data-driven',...
    'Location', 'northeast', 'Interpreter', 'latex')
title('Lyapunov function decrease')
xlabel('t')
ylabel('lyap')
grid on

ax = gca; % Get the current axes handle

ax.Title.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Legend.FontSize = 16;


%% Figure 4: Consistency set shrinkage
if size(X_rh,1) <= 2
    compare_consistency
end