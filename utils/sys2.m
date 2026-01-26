function dxdt = sys2(t,x,sys,N,eps,type)       % 2nd order
global x_glo
global index
u = randn(N,1);
x_glo(1, index) = t;
x_glo(2, index) = x(1);
x_glo(3, index) = x(2);
noise = generate_noise(N,1,eps,type);
dxdt = sys.a*[x(1);x(2)] + sys.b*u + noise;
x_glo(4, index) = dxdt(1);
x_glo(5, index) = dxdt(2);
x_glo(6, index) = u(1);
x_glo(7, index) = u(2);
x_glo(8, index) = noise(1);
x_glo(9, index) = noise(2);
index = index + 1;
end