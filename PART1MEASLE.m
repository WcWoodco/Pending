% %%% MEASLES %%% %
h = 0.5; % day
% Initial conditions
S0 = 990; 
I0 = 10;
R0 = 0; 

a = 0; % day 0 
b = 100; % day 100 (time simulation)
n = (b-a)/h; % step size
B = 2; % transmission rate for measles
g = 0.2; % recovery rate for measles
N = S0 + R0 + I0; % total population (constant)
t = linspace(a,b,n+1);
S = zeros(1,n+1);
I = zeros(1,n+1);
R = zeros(1,n+1);
S(1) = S0;
I(1) = I0;
R(1) = R0;

% Equations
fS = @(S,I) -B*S*I / N; % Susceptible function
fI = @(S,I) B*S*I / N - g*I; % Infected function
fR = @(I) g*I; % Recovered function 

% Fourth order Runge-Kutta method
for i = 1:n
    K1S = fS(S(i),I(i));
    K2S = fS(S(i)+0.5*h, I(i)+0.5*K1S*h);
    K3S = fS(S(i)+0.5*h, I(i)+0.5*K2S*h);
    K4S = fS(S(i+1), I(i)+K3S*h);
    S(i+1) = S(i) + (K1S+2*K2S+2*K3S+K4S)*(h/6);

    K1I = fI(S(i),I(i));
    K2I = fI(S(i) + 0.5*h, I(i) + 0.5*K1I*h);
    K3I = fI(S(i) + 0.5*h, I(i) + 0.5*K2I*h);
    K4I = fI(S(i+1), I(i) + K3I*h);
    I(i+1) = I(i) + (K1I+2*K2I+2*K3I+K4I)* (h/6);

    K1R = fR(I(i));
    K2R = fR(I(i) + 0.5*K1R*h);
    K3R = fR(I(i) + 0.5*K2R*h);
    K4R = fR(I(i) + K3R*h);
    R(i+1) = R(i) + (K1R+2*K2R+2*K3R+K4R)*(h/6);
end

plot(t,S,'b-')
grid on
hold on
plot(t,I,'r-')
plot(t,R,'k-')
legend('Susceptible','Infected','Recovered','Location','east')
title 'Measles SIR Model'
xlabel('Day')
ylabel('Population')
