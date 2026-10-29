function [u, Gamma, P, sol] = MBLQR(A,B,Q,R,UU,x,P_old, opts)
% MBLQR  Model-based input-constrained LQR step (A and B known).
%   [u, Gamma, P, sol] = MBLQR(A, B, Q, R, UU, x, P_old, opts) returns the
%   input u = K*x at the current state x with ||UU*u||_2 <= 1.

n = size(B,1);
m = size(B,2);
p = size(UU,1);

P = sdpvar(n);
M = sdpvar(m,n,'full');
W = sdpvar(m);
Gamma = sdpvar(1);
Beta = sdpvar(1);
Tau = sdpvar(1);

Cons = [
    Gamma - trace(Q*P) - trace(W) >= 0;
    [P-Tau*eye(n), A*P+B*M;(A*P+B*M)', P] >= 0;
    [Beta*eye(p), UU*M; (UU*M)', P] >= 0;
    [Beta, Beta*x'; Beta*x, P] >= 0;
    [W, sqrtm(R)*M; (sqrtm(R)*M)', P] >= 0;
    Tau >= 1;
    P >= P_old;
];

sol = optimize(Cons, Gamma, opts);

P = value(P);
M = value(M);
W = value(W);
K = M/P;
u = K*x;
Gamma = value(Gamma);
Tau = value(Tau);

end