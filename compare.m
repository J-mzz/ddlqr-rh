close all; clear; clc
yalmip('clear')
rng(1)

%% Compare DD-LQR, Input-Constrained DD-LQR, and RH IC DD-LQR 
% Algorithm 2 in ref [2], Algorithm 2, and Algorithm 3, respectively

%% System dynamics and data

% number of initial measurements
T = 5;      % # of samples

% noise level
eps = 1e-1; % noise level
tol = 1e-6; % SDP tolerance: P>0 --> P>=tol*I

% ground truth system
A = [0, -.99; .99, 0];
B = [0; 1];


% generate auxiliary state-space model 
n = size(A, 1);      % # of state
m = size(B, 2);      % # of input

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

clear sysd

%% Hyper-parameters

% relaxation
order = 1;      % P-satz relaxation order
relax = 1;      % relaxation order in Remark 2, ref [2]
partial = 0;    % 1 if know some entries

% yalmip
opts = sdpsettings('solver', 'mosek', 'verbose', 0);

% experiment setting
num_iter = 10;   % # of iterations
tau_len = 50;   % # of Tau for line search, default of logspace is 50
NOISE = eps * (2*rand(n, num_iter) -1);     % fixed random noise

save_flag = false;   % whether data need to be saved 

% initial values 
Tau_init = 5; %1e2;
x_init = 5*(2*rand(n,1) -1); %1*randn(n,1);   % initial point
P_init = zeros(2);

%% Case 1: DD-LQR, Algorithm 2 in ref [2]

% initialization
x_old = x_init;

% log
X_log = x_old;
U_log = [];
Tau_Gamma_log = [];

tic
% DD-LQR
[K,Gamma,sol] = DDLQR_D(train_data,T,Q,R,eps,order,tol,relax,opts,partial);
if sol.problem ~= 0
    warning('solver error with %d', sol.problem)
end

% iteration
for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)
    
    % update dynamic
    u = K*x_old;
    noise = NOISE(:, iter);
    x_old = A*x_old + B*u + noise;

    % check 
    lost = x_old'*x_old + u*u;

    % log 
    X_log = [X_log, x_old];
    U_log = [U_log, u];
end
t_lqr = toc;
loss_lqr = sum(diag(X_log*X_log')) + U_log*U_log'

% save 
if save_flag
    save('ddlqr',"X_log","U_log","Gamma")
end

%% Case 2: Input-constrained DD-LQR, Algorithm 2

% initialization
Tau_old = Tau_init;
x_old = x_init;
P_old = P_init;

% log
X_log = x_old;
U_log = [];
Tau_Gamma_log = [];
Tau_log = [];
P_log = cell(1, num_iter);

tic
% iteration
for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)

    % linear search over Tau
    Tau_Gamma = Inf;
    K = [];

    % use parfor here for large logspace simulation
    TG_result = zeros(1,tau_len);
    K_result = zeros(tau_len, n);
    tau_range = logspace(0, log10(Tau_old), tau_len);
    P_result = cell(1,tau_len);

    parfor i = 1:tau_len 
        Tau_i = tau_range(i);
        [K_temp,Tau_Gamma_temp,sol,P] = ...
            DDLQR_D_refining(train_data,T,Q,R,eps,order,tol,relax,opts,partial, Tau_i, x_old, P_old);
        
        if sol.problem ~= 0
            warning('solver error with %d at tau=%3f', sol.problem, Tau_i)
        end
        
        TG_result(i) = Tau_Gamma_temp;
        K_result(i,:) = K_temp;
        P_result{i} = P;
    end
    
    % find the min tau_gamma 
    [Tau_Gamma, idx] = min(TG_result);
    Tau = tau_range(idx);
    K = K_result(idx,:);
    P = P_result{idx};

    % update dynamic
    u = K*x_old;
    noise = NOISE(:, iter);
    x_old = A*x_old + B*u + noise;
    Tau_old = Tau;
    P_old = P;

    % check 
    lost = x_old'*x_old + u*u;

    % log 
    X_log = [X_log, x_old];
    U_log = [U_log, u];
    Tau_log = [Tau_log, Tau];
    Tau_Gamma_log = [Tau_Gamma_log, Tau_Gamma];
    P_log{iter} = P;
end
t_ic = toc;
loss_ic = sum(diag(X_log*X_log')) + U_log*U_log'

% save 
if save_flag
    save('ddlqr_ic',"Tau_Gamma_log","Tau_log","U_log","X_log","data")
end

%% Case 3: Receding-Horizon Input-Constrained DD-LQR, Algorithm 3

% initialization
Tau_old = Tau_init;
x_old = x_init;
P_old = P_init;

% log
X_log = x_old;
U_log = [];
Tau_Gamma_log = [];
Tau_log = [];
P_log = cell(1, num_iter);

tic
for iter = 1:num_iter
    fprintf('iteration = %d\n', iter)

    % linear search over Tau
    Tau_Gamma = Inf;
    K = [];

    % use parfor here for large logspace and long simulation
    TG_result = zeros(1,tau_len);
    K_result = zeros(tau_len,n);
    tau_range = logspace(0, log10(Tau_old),tau_len);
    P_result = cell(1,tau_len);

    parfor i = 1:tau_len
        Tau_i = tau_range(i);
        [K_temp,Tau_Gamma_temp,sol, P] = ...
            DDLQR_D_refining(data,T,Q,R,eps,order,tol,relax,opts,partial, Tau_i, x_old, P_old);

        if sol.problem ~= 0
            warning('solver error with %d at tau=%3f', sol.problem, Tau_i)
        end

        TG_result(i) = Tau_Gamma_temp;
        K_result(i,:) = K_temp;
        P_result{i} = P;
    end
    
    % find the min tau_gamma
    [Tau_Gamma, idx] = min(TG_result);
    Tau = tau_range(idx);
    K = K_result(idx,:);
    P = P_result{idx};

    % update dynamic
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
    Tau_Gamma_log = [Tau_Gamma_log, Tau_Gamma];
    P_log{iter} = P;
end
t_rh = toc;
loss_rh = sum(diag(X_log*X_log')) + U_log*U_log'

% save
if save_flag
    save('ddlqr_rh',"Tau_Gamma_log","Tau_log","U_log","X_log","data")
end


%% Save the workspace
save('alldata')