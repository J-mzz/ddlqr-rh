function Gamma = H2_continuous(sys, K, type)
N = size(sys.a,1);
switch type
    case 'Slemma'
        E = eye(N);   % matrix of noise
        A_K = sys.a+sys.b*K;
        C_K = sys.c+sys.d*K;
        P = lyap(A_K',C_K'*C_K);
        Gamma = trace(E'*P*E);   %
    case 'Claudio'
        Wx = sys.c'*sys.c;
        Wu = sys.d'*sys.d;
        Phi = sys.a+sys.b*K;
        W = lyap(Phi,eye(N));
        Gamma = trace((Wx+K'*Wu*K)*W);
end 