function [p_out, cons, Gram] = psatz_scalar(X, order, vars, var_n, relax, n, m)
% X         eq/ineq constraints
% order     degree of psatz
% vars      varriables of psatz
% var_n     # of variables
% relax     relaxation in Remark 2
% n         # of state
% m         # of input

L_i = length(X.ineq);       % # of inequality
d_i = degree(X.ineq);       % degree of inequality

Gram = cell(1+L_i, 1);
[pow, ~] = momentPowers(0, var_n, order);   % sigma0
vect_0 = recovermonoms(pow, vars);
% vect_0(1) = [];
if relax == 1
    if var_n == 3*n+m+1
        E1 = sdpvar(1);
        E2 = sdpvar(2*n);
        E3 = sdpvar(n+m);
    else
        E1 = sdpvar(1);
        E2 = sdpvar(n);
        E3 = sdpvar(n+m);
    end
    Gram{1} = blkdiag(E1,E2,E3);
%         cons = (Gram{1} >= 0):'Gram_0';
    cons = [(E1 >= 0):'E1';(E2 >= 0):'E2'; (E3 >= 0):'E3'];
    vect_0(1) = [];
    p_out = vect_0'*Gram{1}*vect_0;
else
    l_0 = length(vect_0);
    Gram_0 = sdpvar(l_0);
    p_out = vect_0'*Gram_0*vect_0;
    cons = (Gram_0 >= 0):'Gram_0';
    Gram{1} = Gram_0;
end

[pow, ~] = momentPowers(0, var_n, floor((2*order-d_i)/2));
vect_i = recovermonoms(pow, vars);
l_i = length(vect_i);
for i = 1:L_i       % sigmai
    Gram{i+1} = sdpvar(l_i);
    p_out = p_out + vect_i'*Gram{i+1}*vect_i*X.ineq(i);
    cons = [cons; (Gram{i+1} >= 0):'Gram_i'];
end

end