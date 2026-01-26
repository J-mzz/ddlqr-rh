function Gamma = compute_H2(sysd,Q,R,type,opts)
A = sysd.a;
B = sysd.b;
n = size(A,1);
m = size(B,2);
tol = 1e-6;
Y = sdpvar(n,n);
M = sdpvar(m,n);
X = sdpvar(m,m);
Gamma = sdpvar(1,1);

if type == 'c'
    C1 = -trace(Q*Y)-trace(X)+Gamma;
    C2 = -Y*A'-A*Y-M'*B'-B*M-eye(n);
    C3 = [X sqrtm(R)*M; M'*sqrtm(R) Y];
elseif type == 'd'
    C1 = -trace(Q*Y)-trace(X)+Gamma;
    C2 = [Y-eye(n) A*Y + B*M; Y*A' + M'*B' Y];
    C3 = [X sqrtm(R)*M; M'*sqrtm(R) Y];
end
LMI = [C1 >= tol; C2 >= tol; C3 >= tol; Y >= tol];

optimize(LMI,Gamma,opts);
Gamma = value(Gamma);
end