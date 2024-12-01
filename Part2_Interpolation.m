% %%% SEASONAL INFLUENZA %%% %
h = 2; % day
% Initial conditions %
S0 = 990; 
I0 = 10;
R0 = 0;

a = 0; % day 0 
b = 100; % day 100 (time simulation)
n = (b-a)/h; % step size
B = 0.3; % transmission rate for seasonal influenza
g = 0.1; % recovery rate for seasonal influenza
N = S0 + R0 + I0; % total population (constant)
t = linspace(a,b,n+1); % x-axis 0 - 100 days
S = zeros(1,n+1);
I = zeros(1,n+1);
R = zeros(1,n+1);
S_1 = zeros(1,2*n+1);
I_1 = zeros(1,2*n+1);
R_1 = zeros(1,2*n+1);
S(1) = S0;
I(1) = I0;
R(1) = R0;
S_1(1) = S0;
I_1(1) = I0;
R_1(1) = R0;

% Equations %
fS = @(S,I) -B*S*I / N; % Susceptible function
fI = @(S,I) B*S*I / N - g*I; % Infected function
fR = @(I) g*I; % Recovered function 

% Fourth Order Runge-Kutta Method for h=2 %
for i = 1:n
    K1S = fS(S(i),I(i));
    K2S = fS(S(i)+0.5*h, I(i)+0.5*K1S*h);
    K3S = fS(S(i)+0.5*h, I(i)+0.5*K2S*h);
    K4S = fS(S(i+1), I(i)+K3S*h);
    S(i+1) =S(i) + (K1S+2*K2S+2*K3S+K4S)*(h/6);

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

% Fourth Order Runge-Kutta Method for h=1%
for i = 1:2*n
    K1S_1 = fS(S_1(i),I_1(i));
    K2S_1 = fS(S_1(i)+0.5*(h/2), I_1(i)+0.5*K1S_1*(h/2));
    K3S_1 = fS(S_1(i)+0.5*(h/2), I_1(i)+0.5*K2S_1*(h/2));
    K4S_1 = fS(S_1(i+1), I_1(i)+K3S_1*(h/2));
    S_1(i+1) = S_1(i) + (K1S_1+2*K2S_1+2*K3S_1+K4S_1)*((h/2)/6);

    K1I_1 = fI(S_1(i),I_1(i));
    K2I_1 = fI(S_1(i) + 0.5*(h/2), I_1(i) + 0.5*K1I_1*(h/2));
    K3I_1 = fI(S_1(i) + 0.5*(h/2), I_1(i) + 0.5*K2I_1*(h/2));
    K4I_1 = fI(S_1(i+1), I_1(i) + K3I_1*(h/2));
    I_1(i+1) = I_1(i) + (K1I_1+2*K2I_1+2*K3I_1+K4I_1)*((h/2)/6);

    K1R_1 = fR(I_1(i));
    K2R_1 = fR(I_1(i) + 0.5*K1R_1*(h/2));
    K3R_1 = fR(I_1(i) + 0.5*K2R_1*(h/2));
    K4R_1 = fR(I_1(i) + K3R_1*(h/2));
    R_1(i+1) = R_1(i) + (K1R_1+2*K2R_1+2*K3R_1+K4R_1)*((h/2)/6);
end

V_int_S_Linear=zeros(1,n+1); %vector for interpolated values
V_model_S_Linear=zeros(1,n+1); %vector for mode values;
V_int_R_Linear=zeros(1,n+1);
V_model_R_Linear=zeros(1,n+1);
V_int_I_Linear=zeros(1,n+1);
V_model_I_Linear=zeros(1,n+1);

% Newton's linear interpolation Method %
for k = 1:n
    V_int_S_Linear(k)=S(k)+((S(k+1)-S(k))/((k+2)-k))*((k+1)-k);
    V_model_S_Linear(k)=S_1(2*k-1);
    V_int_R_Linear(k)=R(k)+((R(k+1)-R(k))/((k+2)-k))*((k+1)-k);
    V_model_R_Linear(k)=R_1(2*k-1);
    V_int_I_Linear(k)=I(k)+((I(k+1)-I(k))/((k+2)-k))*((k+1)-k);
    V_model_I_Linear(k)=I_1(2*k-1);
end

V_int_S_Quadratic=zeros(1,n+1); %vector for interpolated values
V_model_S_Quadratic=zeros(1,n+1); %vector for mode values;
V_int_R_Quadratic=zeros(1,n+1);
V_model_R_Quadratic=zeros(1,n+1);
V_int_I_Quadratic=zeros(1,n+1);
V_model_I_Quadratic=zeros(1,n+1);


% Newton Quadratic Interpolation %
for k = 1:n-1
   b0_S=S(k);
   b1_S=((S(k+1)-S(k))/((k+2)-k));
   b11_S=((S(k+2)-S(k+1))/((k+4)-(k+2)));
   b2_S=(b11_S-b1_S)/(k+4-k);
   V_int_S_Quadratic(k)=b0_S+b1_S*(k+1-k)+b2_S*(k+1-k)*(k+1-k+2);
   V_model_S_Quadratic(k)=S_1(2*k-1);
   b0_I=I(k);
   b1_I=((I(k+1)-I(k))/((k+2)-k));
   b11_I=((I(k+2)-I(k+1))/((k+4)-(k+2)));
   b2_I=(b11_I-b1_I)/(k+4-k);
   V_int_I_Quadratic(k)=b0_I+b1_I*(k+1-k)+b2_I*(k+1-k)*(k+1-k+2);
   V_model_I_Quadratic(k)=I_1(2*k-1);
   b0_R=R(k);
   b1_R=((R(k+1)-R(k))/((k+2)-k));
   b11_R=((R(k+2)-R(k+1))/((k+4)-(k+2)));
   b2_R=(b11_R-b1_R)/(k+4-k);
   V_int_R_Quadratic(k)=b0_R+b1_R*(k+1-k)+b2_R*(k+1-k)*(k+1-k+2);
   V_model_R_Quadratic(k)=R_1(2*k-1);
end

%Error Norm El2 for I(t), R(t), S(t)%
N_tot = n;
Idiff_Linear = V_int_I_Linear-V_model_I_Linear;
Rdiff_Linear = V_int_R_Linear-V_model_R_Linear;
Sdiff_Linear = V_int_S_Linear-V_model_S_Linear;
EL2_I_Linear = sqrt((sum((Idiff_Linear).^2))/(N_tot+1));
EL2_R_Linear = sqrt((sum((Rdiff_Linear).^2))/(N_tot+1));
EL2_S_Linear = sqrt((sum((Sdiff_Linear).^2))/(N_tot+1));
Idiff_Quadratic = V_int_I_Quadratic-V_model_I_Quadratic;
Rdiff_Quadratic = V_int_R_Quadratic-V_model_R_Quadratic;
Sdiff_Quadratic = V_int_S_Quadratic-V_model_S_Quadratic;
EL2_I_Quadratic = sqrt((sum((Idiff_Quadratic).^2))/(N_tot+1));
EL2_R_Quadratic = sqrt((sum((Rdiff_Quadratic).^2))/(N_tot+1));
EL2_S_Quadratic = sqrt((sum((Sdiff_Quadratic).^2))/(N_tot+1));















