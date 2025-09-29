function data = generate_continuous(n1,n2,sys,N,S,eps,tspan,type)
global x_glo
global index
tt = linspace(0,1,S);   % equal distance time stamp
ind = zeros(S,1);
% define data
data.Xm = cell(n1,n2);
data.U = cell(n1,n2);
data.Xp = cell(n1,n2);
data.noise = cell(n1,n2);
for i1= 1:n1
    for i2 = 1:n2
        % data generation
        index = 1;
        initial = randn(N,1);    % initial state
        [t, ~] = ode15s(@(t,x) sys2(t,x,sys,N,eps(i1),type), tspan, initial);
        % Equidistant sampling to get x,u,d
        for k = 1:S
            [~, ind(k)] = min(abs(x_glo(1,:)-tt(k)));    % closest to time stamp
            data.Xm{i1,i2}(:,k) = x_glo(2:3,ind(k));
            data.Xp{i1,i2}(:,k) = x_glo(4:5,ind(k));
            data.U{i1,i2}(:,k) = x_glo(6:7,ind(k));
            data.noise{i1,i2}(:,k) = x_glo(8:9,ind(k));
        end
    end
end
end