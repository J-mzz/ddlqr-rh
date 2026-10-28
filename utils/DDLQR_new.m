function [K,Gamma,sol, P, Tau, M, W] = DDLQR_new(data,T,Q,R,eps,order,tol,relax,opts,partial, Tau_old, x, P_old, UU)
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
% l = 3*n+m+1 + 2*(n+1);
l = (1) + (2*n) + (p+n) + (n+1) + (m+n);
%   zz1    zz2     zz3     zz4     zz5


% define variables
z = sdpvar(l,1);    % indeterminate x in paper
P = sdpvar(n);
M = sdpvar(m,n,'full');
W = sdpvar(m);
Gamma = sdpvar(1);
Beta = sdpvar(1);
Tau = sdpvar(1);

% constraint term   eq.(21)
Z = z*z';     

% selectors
zz1 = [1, zeros(1,l-1)] * Z * [1, zeros(1,l-1)]';                                % 1
zz2 = [zeros(2*n, 1), eye(2*n), zeros(2*n, l-1-2*n)] * Z * ...                   % 2n
      [zeros(2*n, 1), eye(2*n), zeros(2*n, l-1-2*n)]';
zz3 = [zeros(p+n, 1+2*n), eye(p+n), zeros(p+n, (n+1)+(m+n))] * Z * ...           % p+n
      [zeros(p+n, 1+2*n), eye(p+n), zeros(p+n, (n+1)+(m+n))]';
zz4 = [zeros(n+1, 1+2*n+p+n), eye(n+1), zeros(n+1, m+n)] * Z * ...               % n+1
      [zeros(n+1, 1+2*n+p+n), eye(n+1), zeros(n+1, m+n)]';
zz5 = [zeros(m+n, l-(m+n)), eye(n+m)] * Z * [zeros(m+n, l-(m+n)), eye(n+m)]';    % m+n

% cross-term for Xi*Upsilon - 2*PMZ = 0
St = [eye(n), zeros(n)];   % n x 2n (top selector inside 2n block)
Sb = [zeros(n), eye(n)];   % n x 2n (bottom selector inside 2n block)
PMZ = [P; M] * (Sb * zz2 * St');

C1 = (trace(Q*P) + trace(W) - Gamma)            * zz1;
C2 = trace(-[W sqrtm(R)*M; M'*sqrtm(R) P]       * zz5);
C3 = trace(-[P-Tau*eye(n) zeros(n); zeros(n) P] * zz2);
C4 = trace(-[Beta*eye(p), UU*M; M'*UU', P]      * zz3);
C5 = trace(-[Beta, Beta*x'; Beta*x, P]          * zz4);


% psatz of Upsilon variables E_o, E_1
Upsilon = zeros(2*T,n,'like',sdpvar);
Gram = cell(2*n*T,1);
Cons = [];
CON = struct('ineq', [1-z'*z]);
% CON = struct('ineq', []);
for i = 1:2*T    % psatz of each entry of Upsilon
    for j = 1:n
        [Upsilon(i,j), cons_F, Gram{i,1}] = psatz_scalar_refining(CON, order, z, l, relax, n, m, p);
        Cons = [Cons; cons_F];
    end
end

% psatz of whole polynomial   variables F_o, F_1
[p_out, cons_E] = psatz_scalar_refining(CON, order, z, l, relax, n, m, p);
Cons = [Cons; cons_E];
Cons = [Cons; (coefficients(-C1-C2-C3-C4-C5-trace(D*Upsilon) - p_out, z) == 0):'eq'];

% cons for kd == 0, Y > 0, X > 0
kd = coefficients(Xi*Upsilon-2*PMZ,z);
Cons = [Cons; kd == 0; P >= eye(n); W >= tol*eye(m); Gamma >= 0];
Cons = [Cons; P >= P_old];
% Cons = [Cons; 1e6*eye(n) >= P];

% Cons = [Cons; Tau <= Tau_old-tol];
Cons = [Cons; Tau >= 1];
sol = optimize(Cons, Gamma, opts);

% extract solution
P = value(P);
M = value(M);
W = value(W);
K = M/P;
Gamma = value(Gamma);
Tau = value(Tau);

end