function [u_seq, sol] = deepc_solve(Up,Uf,Yp,Yf, ...
        u_ini,y_ini, Q,R, lambda_g,lambda_y, m,p,N, UU, opts)
import yalmip.*
sH  = size(Up,2);
g   = sdpvar(sH,1);
u   = sdpvar(m*N,1);
y   = sdpvar(p*N,1);
s_y = sdpvar(numel(y_ini),1);
s_up = sdpvar(numel(u_ini),1);

Cons = [];                             % 13: should use exact past inputs, 
Cons = [Cons, Up*g == u_ini + s_up];   %     slack for the 0 padding 
Cons = [Cons, Yp*g == y_ini + s_y];    % 14: slack only on past outputs
Cons = [Cons, Uf*g == u, Yf*g == y];

warning off
for k=1:N
    uk = u((k-1)*m+1:k*m);
    Cons = [Cons, norm(UU*uk,2) <= 1];
end
warning on

Qbar = kron(speye(N), Q);
Rbar = kron(speye(N), R);
r    = zeros(p*N,1);
obj  = (y-r)'*Qbar*(y-r) + u'*Rbar*u + lambda_g*(g'*g) + lambda_y*(s_y'*s_y)+ lambda_y*(s_up'*s_up);

sol = optimize(Cons, obj, opts);
u_seq = value(u);
end