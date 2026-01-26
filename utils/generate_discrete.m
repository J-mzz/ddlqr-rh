function data = generate_discrete(sysd,T,eps,type)
% define data
n = size(sysd.a,1);
m = size(sysd.b,2);
x = randn(n,1);
u = randn(m,T);
Noise = generate_noise(n,T,eps,type);
for i = 1:T
    x(:,i+1) = sysd.a*x(:,i) + sysd.b*u(:,i) + Noise(:,i);
end
Xp = x(:, 2:end);
Xm = x(:, 1:end-1);

data.Xm = Xm;
data.U = u;
data.Xp = Xp;
data.Noise = Noise;
end