function noise = generate_noise(n,T,eps,type)
switch type
    case 'gaussian'                     % eps is the standard deviation
        noise = normrnd(0,eps,[n T]);
    case 'inf_norm'
        noise = eps*(2*rand(n,T) - 1);  % eps is the upper bound
    case 'constant'
        noise = randi(2,n,T)-1;
        noise(~noise) = -1;
        noise = eps*noise;
end