function plot_results(folder)
% PLOT_RESULTS  Compare DeePC, DD-LQR, and the input-constrained DD-LQR.
%   plot_results(folder) loads the results saved by ddlqr_2d.m or
%   ddlqr_4d.m in folder and plots the trajectories, the control inputs,
%   the worst-case upper bound, the Lyapunov function, and (2D only) the
%   consistency sets with and without receding horizon.
%
%   plot_results('example_2d') and plot_results('example_4d') plot the
%   results reported in the paper.
%
%   For more than one input, the plotted input is ||UU*u||_2.

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'utils')))
if ~isfolder(folder)
    folder = fullfile(root, folder);   % allow paths relative to the repo
end

close all

S = load(fullfile(folder, 'alldata'), 'UU');
UU = S.UU;

%% load data
S = load(fullfile(folder, 'mdlqr_ic'));
X_md = S.X_log;
U_md = input_norm(UU, S.U_log);
G_md = S.Gamma_log;
V_md = S.V_log;

t = length(U_md);

S = load(fullfile(folder, 'deepc'));
X_dpc = S.X_log;
U_dpc = input_norm(UU, S.U_log);

S = load(fullfile(folder, 'ddlqr'));
X_lqr = S.X_log;
U_lqr = input_norm(UU, S.U_log);
G_lqr = S.Gamma * ones(1, t);
V_lqr = S.V_log;

S = load(fullfile(folder, 'ddlqr_ic'));
X_ic  = S.X_log;
U_ic  = input_norm(UU, S.U_log);
T_ic  = S.Tau_log;
G_ic  = S.Gamma_log;
V_ic  = S.V_log;

S = load(fullfile(folder, 'ddlqr_rh'));
X_rh  = S.X_log;
U_rh  = input_norm(UU, S.U_log);
T_rh  = S.Tau_log;
G_rh  = S.Gamma_log;
V_rh  = S.V_log;


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


%% Figure 1-1: Trajectories, coordinates 3 and 4 (4D example)

if size(X_rh, 1) > 2

    b = 3;

    figure()
    plot(X_lqr(3,:), X_lqr(4,:), '*--', 'LineWidth', 1, 'Color', color1)
    hold on
    plot(X_md(3,:),  X_md(4,:),  '*--', 'LineWidth', 1, 'Color', color3)
    plot(X_rh(3,:),  X_rh(4,:),  '*--', 'LineWidth', 1, 'Color', color2)
    plot(X_dpc(3,:), X_dpc(4,:), '*--', 'LineWidth', 1, 'Color', color4)
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

%% Figure 2: Control

figure()
plot(1:t, U_lqr,'*--', 'LineWidth', 1, 'Color', color1)
hold on
plot(1:t, U_md, '*--', 'LineWidth', 1, 'Color', color3)
plot(1:t, U_rh, '*--', 'LineWidth', 1, 'Color', color2)
plot(1:t, U_dpc,'*--', 'LineWidth', 1, 'Color', color4)
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


%% Figure 4: Consistency set shrinkage (2D example)
if size(X_rh,1) <= 2
    plot_consistency_sets(folder)
end

end


function u = input_norm(UU, U_log)
% ||UU*u_k||_2 for each step if there is more than one input
if size(UU,1) > 1
    u = vecnorm(UU*U_log, 2);
else
    u = U_log;
end
end
