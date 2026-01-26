function [K,Gamma,sol,P,Tau, M, W] = DDLQR_scherer(data,T,Q,R,eps,order,tol,type,opts, x, P_old, UU)
%% Solve using Scherer P-satz
% extract data
Xm = data.Xm(:,1:T);
U = data.U(:,1:T);
Xp = data.Xp(:,1:T);

n = size(Xm,1);   % size of state
m = size(U,1);    % size of input
p = size(UU,1);

% define variables
A = sdpvar(n,n,'full');
B = sdpvar(n,m,'full');
P = sdpvar(n);
M = sdpvar(m,n,'full');
W = sdpvar(m);
eta = sdpvar(1);
Gamma = sdpvar(1);
Tau = sdpvar(1);
Beta = sdpvar(1);

C1 = -trace(Q*P)-trace(W)+Gamma;
C2 = [P-Tau*eye(n), A*P + B*M; P*A' + M'*B', P];
C3 = [W sqrtm(R)*M; M'*sqrtm(R) P];
C4 = [Beta*eye(p), UU*M; M'*UU', P];
C5 = [Beta, Beta*x'; Beta*x, P];

% F = blkdiag(C1,C2,C3);
% F = F-eta*eye(size(F));
F = C2-eta*eye(size(C2));      % only C2 involves A,B
gp = eps - (Xp - A*Xm - B*U);
gm = eps + (Xp - A*Xm - B*U);

Cons = struct('ineq', [gp(:);gm(:)], 'eq', []);
vars = [A(:);B(:)];
var_count = length(vars);
[~, cons] = psatz_matrix_refining(F, Cons, order, vars, var_count);
cons = [cons; P >= eye(n); eta >= tol; W >= tol*eye(m)];
cons = [cons; C1 >= tol; C3 >= tol; C4 >= tol; C5 >= tol];
cons = [cons; P >= P_old; Tau >= 1];

sol = optimize(cons,Gamma,opts);

P = value(P);
M = value(M);
W = value(W);
K = M/P;
Gamma = value(Gamma);
Tau = value(Tau);

end