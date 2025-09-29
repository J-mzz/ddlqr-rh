function Gamma = H2_discrete(sysd, K, type)
N = size(sysd.a,1);
switch type
    case 'Slemma'
        E = eye(N);   % matrix of noise
        A_K = sysd.a+sysd.b*K;
        C_K = sysd.c+sysd.d*K;
        P = dlyap(A_K',C_K'*C_K);
        Gamma = trace(E'*P*E);   %
    case 'Claudio'
        Wx = sysd.c'*sysd.c;
        Wu = sysd.d'*sysd.d;
        Phi = sysd.a+sysd.b*K;
        W = dlyap(Phi,eye(N));
        Gamma = trace((Wx+K'*Wu*K)*W);
end 