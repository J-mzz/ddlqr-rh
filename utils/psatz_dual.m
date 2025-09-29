function [p_out, cons, Gram, Mu] = psatz_dual(p, X, order, vars, var_n, n, T, type)
d_i = degree(X.ineq);       % degree of inequality
d_e = degree(X.eq);         % degree of equality

Gram = cell(1+2*n*T, 1);
Mu = cell(n*(T-1), 1);

[pow, ~] = momentPowers(0, var_n, order);   % sigma0
vect_0 = recovermonoms(pow, vars);
l_0 = length(vect_0);
Gram_0 = sdpvar(l_0);
p_out = vect_0'*Gram_0*vect_0;
cons = (Gram_0 >= 0):'Gram_0';
Gram{1} = Gram_0;

if mod(d_i,2) == 0
    [pow, ~] = momentPowers(0, var_n, floor((2*order-d_i)/2));
else
    [pow, ~] = momentPowers(0, var_n, floor((2*order-d_i-1)/2));
end
vect_zeta = recovermonoms(pow, vars);
l_zeta = length(vect_zeta);
for i = 1:n*T     % inequality positive
    Gram{i+1} = sdpvar(l_zeta);
    cons = [cons; (Gram{i+1} >= 0):'Gram_zeta_p'];
    zeta_p(i,1) = vect_zeta'*Gram{i+1}*vect_zeta;
    p_out = p_out + zeta_p(i,1)*X.ineq(i);
end
zeta_p = reshape(zeta_p, n, T);

for i = 1:n*T     % inequality negative
    Gram{i+1+n*T} = sdpvar(l_zeta);
    cons = [cons; (Gram{i+1+n*T} >= 0):'Gram_zeta_n'];
    zeta_n(i,1) = vect_zeta'*Gram{i+1+n*T}*vect_zeta;
    p_out = p_out + zeta_n(i,1)*X.ineq(i+n*T);
end
zeta_n = reshape(zeta_n, n, T);

[pow, ~] = momentPowers(0, var_n, 2*order-d_e);
vect_mu = recovermonoms(pow, vars);
l_mu = length(vect_mu);
for i = 1:n*(T-1)      % mu
    Mu{i} = sdpvar(l_mu,1);
    mu(i,1) = vect_mu'*Mu{i};
    p_out = p_out + mu(i,1)*X.eq(i);
end
mu = reshape(mu, n, T-1);

switch type
    case 'no_prior'
        A = reshape(vars(1:n^2),n,n);
        B = vars(1+n^2:end);
    case 'prior'
        A = [vars(1) 0.3968; vars(2) 1.0388];
        B = [0.4192; 0.6852];
    case 'prior2'
        A = [vars(1) 0.3968; vars(2) 1.0388];
        B = [vars(3:4)];
end
cons = [cons; (coefficients(p - p_out, vars) == 0):'eq'];
cons = [cons; (coefficients(zeta_p(:,1) - zeta_n(:,1) - A'*mu(:,1), vars) == 0):'eq_extra'];
for i = 2:T-1
    cons = [cons; (coefficients(zeta_p(:,i) - zeta_n(:,i) - (A'*mu(:,i) - mu(:,i-1)), vars) == 0):'eq_extra'];
end
cons = [cons; (coefficients(zeta_p(:,end) - zeta_n(:,end) - (-mu(:,end)), vars) == 0):'eq_extra'];
end