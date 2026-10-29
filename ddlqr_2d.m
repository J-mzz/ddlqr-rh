%% 2D example of "Robust Data-Driven Receding-Horizon Control for LQR With Input Constraints"
%
% System (only used to generate data and simulate):
%   x_{k+1} = A*x_k + B*u_k + w_k,  A = [0 -0.99; 0.99 0],  B = I,
%   ||w_k||_inf <= eps,  input constraint ||UU*u_k||_2 <= 1
%
% Compares, over num_iter closed-loop steps from x_init:
%   - model-based input-constrained LQR (Algorithm 1, knows A and B)
%   - DeePC [11]
%   - DD-LQR [10]: data-driven LQR without input constraints, designed once
%     from the initial data
%   - input-constrained DD-LQR, consistency set fixed
%   - receding-horizon input-constrained DD-LQR (Algorithm 2), consistency
%     set refined with every new measurement
% Case 4 solves the receding-horizon problem with Scherer's P-satz instead
% of Putinar's, for the comparison in the paper. It is off by default; set
% run_scherer = true to run it.
%
% The results are saved to results/2d and plotted by plot_results.m.
% The results reported in the paper are in example_2d/; plot them with
% plot_results('example_2d').

close all; clear; clc
yalmip('clear')
rng(1)

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'utils')))
out_dir = fullfile(root, 'results', '2d');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% System dynamics and data

% number of initial measurements
T = 5;      % # of samples
T_init = T;

% noise level
eps = 1e-1; % noise level
tol = 1e-6; % SDP tolerance: P>0 --> P>=tol*I

% ground truth system
A = [0, -.99; .99, 0];
B = eye(2);

% control constraints: ||Uu|| <= 1
UU = 0.5*eye(2);

% generate auxiliary state-space model
n = size(A, 1);      % # of state
m = size(B, 2);      % # of input
p = size(UU,1);      % # of output


C = eye(n+m,n);
D = [zeros(n,m); eye(m)];
sysd = ss(A,B,C,D,-1);

% LQR
Q = sysd.C'*sysd.C;     % here we select Q = I, R = I
R = sysd.D'*sysd.D;

% generate initial measurements
train_data = generate_discrete(sysd,T,eps,'inf_norm'); % fixed for record
data = train_data;
data.eps = eps;
train_data.eps = eps;

clear sysd

%% Hyper-parameters

% relaxation
order = 1;      % P-satz relaxation order
relax = 1;      % 1: block-diagonal Gram matrix for sigma_0, 0: full
partial = 0;    % 1 if know some entries

% yalmip
opts = sdpsettings('solver', 'mosek', 'verbose', 0);

% experiment setting
num_iter = 10;   % # of iterations
NOISE = eps*(2*rand(n, num_iter) -1);     % fixed random noise

save_flag = true;   % whether data need to be saved
run_scherer = false;    % also run Case 4 (Scherer's P-satz)

% initial values
Tau_init = 1e2;
x_init = 5*[1;1];
P_init = eps*eps/3 * eye(2);

%% Model-based input-constrained LQR (Algorithm 1)
% initialization
x_old = x_init;
P_old = P_init;

% log
X_log = x_old;
U_log = [];
Gamma_log = [];
P_log = cell(1, num_iter);
T_md = [];
V_log = [];

% iteration
for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)

    [u, Gamma, P, sol] = MBLQR(A,B,Q,R,UU,x_old,P_old,opts);
    T_md = [T_md, sol.solvertime];

    if sol.problem ~= 0
        warning('solver error with %d', sol.problem)
    end

    % update dynamic
    V_log = [V_log, x_old'*inv(P)*x_old];
    noise = NOISE(:, iter);
    x_old = A*x_old + B*u + noise;
    P_old = P;

    % log
    X_log = [X_log, x_old];
    U_log = [U_log, u];
    Gamma_log = [Gamma_log, Gamma];
    P_log{iter} = P;

end
loss_md = trace(Q * (X_log*X_log.')) + trace(R * (U_log*U_log.'))

% save
if save_flag
    save(fullfile(out_dir, 'mdlqr_ic'),"Gamma_log","U_log","X_log","V_log")
end


%% Case 0: DeePC [11]

% initialization
x_old = x_init;

% logs
X_log = x_old;
U_log = [];
T_deepc = [];

% parameters
Tini = 2; Tpre = 2; L = Tini + Tpre;
lambda_g  = 1e-3;
lambda_y  = 1e2;
lambda_ini = 1e-2;

% data buffer
u_buf = data.U;
y_buf = data.Xm;

% initial Hankel
[Hu, Hy] = hankel_block(u_buf, y_buf, L);
Up = Hu(1:m*Tini,:);          Uf = Hu(m*Tini+1:end,:);
Yp = Hy(1:p*Tini,:);          Yf = Hy(p*Tini+1:end,:);

for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)

    % --- build past window (zero padded before enough data) ---
    if size(U_log,2) < Tini
        need_u = Tini - size(U_log,2);
        u_hist_m = [zeros(m, need_u), U_log];
    else
        u_hist_m = U_log(:, end-Tini+1:end);
    end

    if size(X_log,2) <= Tini
        need_y = Tini - (size(X_log,2)-1);
        y_hist_m = [zeros(p, need_y), X_log(:, 1:end-1)];
    else
        y_hist_m = X_log(:, end-Tini:end-1);
    end

    % vectorize for DeePC
    u_ini = reshape(u_hist_m, [], 1);
    y_ini = reshape(y_hist_m, [], 1);

    % --- call standard DeePC solver ---
    [u_seq, sol] = deepc_solve( ...
        Up,Uf,Yp,Yf, u_ini,y_ini, Q,R,lambda_g,lambda_y,m,p,Tpre,UU,opts);
    T_deepc = [T_deepc, sol.solvertime];

    if sol.problem ~= 0
        warning('Solver issue at iter %d: %s', iter, sol.info);
    end

    % --- apply first input & step dynamics ---
    u = u_seq(1:m);
    noise = NOISE(:,iter);
    x_old = A*x_old + B*u + noise;

    % log & update data buffer
    U_log = [U_log, u];
    X_log = [X_log, x_old];
    u_buf = [u_buf, u];
    y_buf = [y_buf, x_old];

    % rebuild Hankels for next iteration
    if size(u_buf,2) >= L
        [Hu, Hy] = hankel_block(u_buf, y_buf, L);
        Up = Hu(1:m*Tini,:); Uf = Hu(m*Tini+1:end,:);
        Yp = Hy(1:p*Tini,:); Yf = Hy(p*Tini+1:end,:);
    end
end
loss_deepc = trace(Q * (X_log*X_log.')) + trace(R * (U_log*U_log.'))

% save
if save_flag
    save(fullfile(out_dir, 'deepc'),"X_log","U_log")
end

%% Case 1: DD-LQR [10]

% initialization
x_old = x_init;

% log
X_log = x_old;
U_log = [];
V_log = [];

% DD-LQR
[K,Gamma,sol,P] = DDLQR_D(train_data,T,Q,R,eps,order,tol,relax,opts,partial);
T_lqr = sol.solvertime;

if sol.problem ~= 0
    warning('solver error with %d', sol.problem)
end

% iteration
for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)

    % update dynamic
    V_log = [V_log; x_old'*inv(P)*x_old];
    u = K*x_old;
    noise = NOISE(:, iter);
    x_old = A*x_old + B*u + noise;

    % log
    X_log = [X_log, x_old];
    U_log = [U_log, u];
end

loss_lqr = trace(Q * (X_log*X_log.')) + trace(R * (U_log*U_log.'))

% save
if save_flag
    save(fullfile(out_dir, 'ddlqr'),"X_log","U_log","Gamma","V_log")
end

%% Case 2: Input-constrained DD-LQR (consistency set fixed)

% initialization
Tau_old = Tau_init;
x_old = x_init;
P_old = P_init;

% log
X_log = x_old;
U_log = [];
Gamma_log = [];
Tau_log = [];
P_log = cell(1, num_iter);
T_ic = [];
V_log = [];

% iteration
for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)

    [K,Gamma,sol,P, Tau, M, W] = ...
        DDLQR_IC(train_data,T,Q,R,eps,order,tol,relax,opts,partial, Tau_old, x_old, P_old, UU);
    T_ic = [T_ic, sol.solvertime];

    if sol.problem ~= 0
        warning('solver error with %d', sol.problem)
    end

    % update dynamic
    V_log = [V_log, x_old'*inv(P)*x_old];
    u = K*x_old;
    noise = NOISE(:, iter);
    x_old = A*x_old + B*u + noise;
    Tau_old = Tau;
    P_old = P;

    % log
    X_log = [X_log, x_old];
    U_log = [U_log, u];
    Tau_log = [Tau_log, Tau];
    Gamma_log = [Gamma_log, Gamma];
    P_log{iter} = P;

end
loss_ic = trace(Q * (X_log*X_log.')) + trace(R * (U_log*U_log.'))

% save
if save_flag
    save(fullfile(out_dir, 'ddlqr_ic'),"Gamma_log","Tau_log","U_log","X_log","data","V_log")
end

%% Case 3: Receding-horizon input-constrained DD-LQR (Algorithm 2)

% initialization
Tau_old = Tau_init;
x_old = x_init;
P_old = P_init;

% log
X_log = x_old;
U_log = [];
Gamma_log = [];
Tau_log = [];
P_log = cell(1, num_iter);
T_rh = [];
V_log = [];

for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)

    [K,Gamma,sol,P, Tau, M, W] = ...
        DDLQR_IC(data,T,Q,R,eps,order,tol,relax,opts,partial, Tau_old, x_old, P_old, UU);
    T_rh = [T_rh, sol.solvertime];

    if sol.problem ~= 0
        warning('solver error with %d', sol.problem)
    end

    % update dynamic
    V_log = [V_log, x_old'*inv(P)*x_old];
    u = K*x_old;
    Tau_old = Tau;
    P_old = P;

    % add the new measurement to the data
    data.Xm = [data.Xm, x_old];
    data.U = [data.U, u];

    noise = NOISE(:, iter);
    x_old = A*x_old + B*u + noise;

    data.Xp = [data.Xp, x_old];
    data.Noise = [data.Noise, noise];
    T = T+1;

    % log
    X_log = [X_log, x_old];
    U_log = [U_log, u];
    Tau_log = [Tau_log, Tau];
    Gamma_log = [Gamma_log, Gamma];
    P_log{iter} = P;

end
loss_rh = trace(Q * (X_log*X_log.')) + trace(R * (U_log*U_log.'))

% save
if save_flag
    save(fullfile(out_dir, 'ddlqr_rh'),"Gamma_log","Tau_log","U_log","X_log","data","V_log")
end

%% Case 4: Receding-horizon input-constrained DD-LQR with Scherer's P-satz
if run_scherer
    T = T_init;
    data = train_data;

    % initialization
    Tau_old = Tau_init;
    x_old = x_init;
    P_old = P_init;

    % log
    X_log = x_old;
    U_log = [];
    Gamma_log = [];
    Tau_log = [];
    P_log = cell(1, num_iter);
    T_sch = [];
    V_log = [];

    for iter = 1:num_iter
        fprintf('iteration = %d\n', iter)

        [K,Gamma,sol,P, Tau, M, W] = ...
            DDLQR_scherer(data,T,Q,R,eps,order,tol,'d',opts, x_old, P_old, UU);
        T_sch = [T_sch, sol.solvertime];

        if sol.problem ~= 0
            warning('solver error with %d', sol.problem)
        end

        % update dynamic
        V_log = [V_log, x_old'*inv(P)*x_old];
        u = K*x_old;
        Tau_old = Tau;
        P_old = P;

        data.Xm = [data.Xm, x_old];
        data.U = [data.U, u];

        noise = NOISE(:, iter);
        x_old = A*x_old + B*u + noise;

        data.Xp = [data.Xp, x_old];
        data.Noise = [data.Noise, noise];
        T = T+1;

        % log 
        X_log = [X_log, x_old];
        U_log = [U_log, u];
        Tau_log = [Tau_log, Tau];
        Gamma_log = [Gamma_log, Gamma];
        P_log{iter} = P;

    end
    loss_sch = trace(Q * (X_log*X_log.')) + trace(R * (U_log*U_log.'))

    % save
    if save_flag
        save(fullfile(out_dir, 'ddlqr_sch'),"Gamma_log","Tau_log","U_log","X_log","data","V_log")
    end
end

%% Save the workspace
save(fullfile(out_dir, 'alldata'))

%% Plot
plot_results(out_dir)


%% Auxiliary functions
function [Hu, Hy] = hankel_block(u, y, L)
[m, T] = size(u);
[p, ~] = size(y);

s = T - L + 1;
Hu = zeros(m*L, s);
Hy = zeros(p*L, s);
for i = 1:s
    Hu(:,i) = reshape(u(:, i:i+L-1), m*L, 1);
    Hy(:,i) = reshape(y(:, i:i+L-1), p*L, 1);
end
end
