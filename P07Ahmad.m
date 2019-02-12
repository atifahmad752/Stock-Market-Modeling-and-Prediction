function [S,mu,sigma] = P07Ahmad(S0,mu0,sigma0,Sday)
i = 1;
N = 390;
dt = (ceil(max(Sday))-floor(min(Sday)))/N;
%for i = 1:20
t(1) = 0;
for l = 2:N+1
    t(l) = t(l-1) + dt;
end

for i = 1:77
w(1) = 0;
sigma1 = sqrt(dt);
for k = 2:N+1
    a(k) = sigma1*randn(1,1);
    w(k) = w(k-1) + a(k);
end


q(1) = S0;
q(1)
for j = 2:6
    q(j) = S0*exp((mu0-(((sigma0).^2)/2)*t(j))+sigma0*w(j));
end

y = q(2:6);


S(5*(i-1)+1:5*(i-1)+5) = y;

S0 = Sday(5*(i-1)+5);

%t = t(2:N+1);
%mu0 = S0*exp(mu0*t(6*i));
%average = S0*exp(mu0*(t(5*((i+1)-1)+5)-t(5*((i+1)-1)+1)));
c = mean(Sday((5*(i-1)+1):(5*(i-1)+5)));

%mu0 = (log(S0/Sday(5*(i-1)+1)))/((5*(i-1)+5)-(5*(i-1)+1));
mu0 = (log(c/Sday(5*(i-1)+1)))/5*dt;

%mu0 = S0*exp(mu0*5*i*dt);
%mu0 = S0*exp(mu0*5*dt);
%average = S0*exp(mu0*5*i*dt);
%mu0 = (average-mu0)/mu0;
mu(i) = mu0;
%sigma0 = S0*exp(mu0*t(6*i))*sqrt(exp((sigma0.^2)*t(6*i))-1);
%sigma0 = S0*exp(mu0*(t(5*((i+1)-1)+5)-t(5*((i+1)-1)+1)))*sqrt(exp((sigma0.^2)*(t(5*((i+1)-1)+5)-t(5*((i+1)-1)+1)))-1);

%v = log(var(Sday(5*(i-1)+1:5*(i-1)+5))*((Sday(5*(i-1)+1)).^(-2))*exp(-2*mu0*(5*dt))+1)/(5*dt);
%v = (log((var(Sday(5*(i-1)+1:5*(i-1)+5))/(Sday(5*(i-1)+1)).^2)+exp(2*mu0*5*dt)))/(5*dt)
%Correct code below
v = (log((var(Sday(5*(i-1)+1:5*(i-1)+5)))/(((Sday(5*(i-1)+1)).^2)*(exp(2*mu0*5*dt)))+1))/(5*dt);
sigma0 = sqrt(v);

%sigma0 = sqrt( (log( (std(Sday(5*(i-1)+1:5*(i-1)+5))/S0*exp(mu0*((5*(i-1)+5)-(5*(i-1)+1)))).^2 + 1)/((5*(i-1)+5)-(5*(i-1)+1))));
%sigma0 = sqrt( (log( (std(Sday(5*(i-1)+1:5*(i-1)+5))/(Sday(5*(i-1)+1))*exp(mu0*((5*(i-1)+5)-(5*(i-1)+1)))).^2 + 1)/(t(5*(i-1)+5)-t(5*(i-1)+1))));

%sigma0 = S0*exp(mu0*5*i*dt)*sqrt(exp((sigma0.^2)*5*i*dt)-1);
%sigma0 = S0*exp(mu0*5*dt)*sqrt(exp((sigma0.^2)*5*dt)-1);
%average2 = S0*exp(mu0*5*i*dt)*sqrt(exp((sigma0.^2)*5*i*dt)-1);
%sigma0 = (average2 - sigma0)/sigma0;
sigma(i) = sigma0;
%S0 = Sday(5*(i-1)+5);
end

t = t(2:386);
Sday = Sday(1:385);
plot(t,Sday,'red');
hold on;
plot(t,S,'blue');
end