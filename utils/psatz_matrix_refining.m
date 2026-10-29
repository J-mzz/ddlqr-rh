function [p_out, cons, Gram] = psatz_matrix_refining(p, X, order, vars, var_n)
% PSATZ_MATRIX_REFINING  Scherer's matrix-valued P-satz certificate that
%   the polynomial matrix p is positive semidefinite on the set given by
%   the constraints in X (used by DDLQR_SCHERER).

l_p = length(p);    % p in paper
L_i = length(X.ineq);
L_mu = length(X.eq);
d_i = degree(X.ineq);

Gram = cell(1+L_i, 1);
Mu = cell(L_mu, 1);

[pow, ~] = momentPowers(0, var_n, order);   % sigma0
vect_0 = recovermonoms(pow, vars);
l_0 = length(vect_0);     % q in paper
temp = kron(vect_0, eye(l_p));
Gram_0 = sdpvar(l_0*l_p);
p_out = temp'*Gram_0*temp;
cons = Gram_0 >= 0;
Gram{1} = Gram_0;

for i = 1:L_i    % sigmai
    if mod(d_i,2) == 0
        [pow, ~] = momentPowers(0, var_n, floor((2*order-d_i)/2));
    else
        [pow, ~] = momentPowers(0, var_n, floor((2*order-d_i-1)/2));
    end
    vect_i = recovermonoms(pow, vars);
    l_i = length(vect_i);
    temp = kron(vect_i, eye(l_p));
    Gram{i+1} = sdpvar(l_i*l_p);
    p_out = p_out + temp'*Gram{i+1}*temp*X.ineq(i);
    cons = [cons; Gram{i+1} >= 0];
end

for i = 1:L_mu    %% not tested yet
    [pow, ~] = momentPowers(0, var_n, 2*order-1);
    vect_mu = recovermonoms(pow, vars);
    l_mu = length(vect_mu);
    Mu{i} = sdpvar(l_mu,1);
    p_out = p_out + vect_mu'*Mu{i}*X.eq(i);
end

cons = [cons; coefficients(p - p_out, vars) == 0];
end

