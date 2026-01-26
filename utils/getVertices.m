function vert = getVertices(data)
% get vertices from data = (xk+1, xk, uk)

Xp = data.Xp;
Xm = data.Xm;
U = data.U;
n = size(Xp,2);
eps = data.eps;

aa = []; bb = [];
for i = 1:n
    aa = [aa; kron(eye(2),Xm(:, i)')];
    bb = [bb; kron(eye(2),U(:, i)')];
end
xi = Xp(:);

N = [aa bb; -aa -bb];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end

e = [eps*ones(size(xi,1),1)+xi; eps*ones(size(xi,1),1)-xi];
[N, e] = nontrivial_constraints(N, e);
vert = sortrows(con2vert(N,e));

end