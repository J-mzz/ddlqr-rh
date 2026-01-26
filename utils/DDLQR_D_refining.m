function [K,Tau_Gamma,sol, P] = DDLQR_D_refining(data,T,Q,R,eps,order,tol,relax,opts,partial, Tau, x, P_old, UU)
% extract data     eq.(17)
Xm = data.Xm(:,1:T);
U = data.U(:,1:T);
Xp = data.Xp(:,1:T);
Xi = [Xm -Xm; U -U];
D = [eps+Xp, eps-Xp];

% size of problem
n = size(Xm,1);   % size of state
m = size(U,1);    % size of input
p = size(UU, 1);  % size of input constraint
l = 3*n+m+1 + 2*(n+1);

% partial info
if partial == 1
    T = T + n;
    L = 30;   % large number as trivial bound
    Xm = [Xm eye(n)];
    U = [U zeros(m,n)];
%     A_partial = [ -0.1660   -0.8153   -0.1616
%         0.4406   -0.6275    0.3704
%         -0.9998   -0.3089   -0.5911];
%     A_plus = [A_partial L*ones(3,2); L*ones(2,5)];
%     A_minus = [-A_partial L*ones(3,2); L*ones(2,5)];
    A_partial = [zeros(1,n-1); eye(n-1)];
    A_plus = [A_partial L*ones(n,1)];
    A_minus = [-A_partial L*ones(n,1)];
    Xi = [Xm -Xm; U -U];
    D = [eps+Xp, A_plus, eps-Xp, A_minus];
end

% define variables
z = sdpvar(l,1);    % indeterminate x in paper
P = sdpvar(n);
M = sdpvar(m,n,'full');
W = sdpvar(m);
Gamma = sdpvar(1);

Beta = sdpvar(1);
% Tau = sdpvar(1);

% Tau is viewed as a constant in each step of linear search

% constraint term   eq.(21)
Z = z*z';     
PMZ = [P;M]*[zeros(n,1+n) eye(n) zeros(n,n+m) zeros(n, 2*(n+1))]*Z*[zeros(n,1) eye(n) zeros(n,n+m+n) zeros(n, 2*(n+1))]';
C1 = (trace(Q*P) + trace(W) - Gamma)*[1 zeros(1,l-1)]*Z*[1 zeros(1,l-1)]';
C2 = trace(-[W sqrtm(R)*M; M'*sqrtm(R) P]*[zeros(n+m,1+2*n) eye(n+m) zeros(n+m, 2*(n+1))]*Z*[zeros(n+m,1+2*n) eye(n+m) zeros(n+m, 2*(n+1))]');
C3 = trace(-[P-1/Tau*eye(n) zeros(n); zeros(n) P]*[zeros(2*n,1) eye(2*n) zeros(2*n,n+m) zeros(2*n, 2*(n+1))]*Z*[zeros(2*n,1) eye(2*n) zeros(2*n,n+m) zeros(2*n, 2*(n+1))]');


C4 = trace(-[Beta*eye(p), UU*M; M'*UU', P]*[zeros(p+n, l-(p+n)-(1+n)), eye(p+n), zeros(p+n, 1+n)]*Z*[zeros(p+n, l-(p+n)-(1+n)), eye(p+n), zeros(p+n, 1+n)]');
C5 = trace(-[Beta, Beta*x'; Beta*x, P]*[zeros(1+n, l-(1+n)), eye(1+n)]*Z*[zeros(1+n, l-(1+n)), eye(1+n)]');

% psatz of Upsilon variables E_o, E_1
Upsilon = zeros(2*T,n,'like',sdpvar);
Gram = cell(2*n*T,1);
Cons = [];
CON = struct('ineq', [1-z'*z]);
% CON = struct('ineq', []);
for i = 1:2*T    % psatz of each entry of Upsilon
    for j = 1:n
        [Upsilon(i,j), cons_F, Gram{i,1}] = psatz_scalar_refining(CON, order, z, l, relax, n, m);
        Cons = [Cons; cons_F];
    end
end

% psatz of whole polynomial   variables F_o, F_1
[p_out, cons_E] = psatz_scalar_refining(CON, order, z, l, relax, n, m);
Cons = [Cons; cons_E];
Cons = [Cons; (coefficients(-C1-C2-C3-C4-C5-trace(D*Upsilon) - p_out, z) == 0):'eq'];

% cons for kd == 0, Y > 0, X > 0
kd = coefficients(Xi*Upsilon-2*PMZ,z);
Cons = [Cons; kd == 0; P >= tol*eye(n); W >= tol*eye(m)];
Cons = [Cons; P >= P_old];
sol = optimize(Cons, Tau*Gamma, opts);

% extract solution
P = value(P);
M = value(M);
K = M/P;
Gamma = value(Gamma);
Tau_Gamma = Tau*Gamma;

end