function [S] = P06Ahmad(S0,mu,sigma)
%Simulate Geometric Brownian Motion on the time interval [0,1] with [0,1]
%equally divided into N = 390 subintervals.
N = 390;
dt = 1/N;
t(1) = 0;
for i = 2:N+1
    t(i) = t(i-1) + dt;
end

w(1) = 0;
sigma1 = sqrt(dt);
for i = 2:N+1
    a(i) = sigma1*randn(1,1);
    w(i) = w(i-1) + a(i);
end

S(1) = S0;
for i = 2:N+1
    S(i) = S0*exp((mu-(((sigma).^2)/2)*t(i))+sigma*w(i));
end

t = t(2:N+1);
S = S(2:N+1);
plot(t,S);
end