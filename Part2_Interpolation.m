% %%% SEASONAL INFLUENZA %%% %
h = 1; % day
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
S_1 = zeros(1,n+1);
I_1 = zeros(1,n+1);
R_1 = zeros(1,n+1);
S(1) = S0;
I(1) = I0;
R(1) = R0;

% Equations %
fS = @(S,I) -B*S*I / N; % Susceptible function
fI = @(S,I) B*S*I / N - g*I; % Infected function
fR = @(I) g*I; % Recovered function 

%% Fourth Order Runge-Kutta Method for h=1 %%
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

S_linear=S;
R_linear=R;
I_linear=I;

%% Making vector into stepsize h=2 %%

for i=1:(n/2)
    S_linear(2*i)=0;
    R_linear(2*i)=0;
    I_linear(2*i)=0;

end

%%Setting up vectors for linear interpolation%%

V_int_S_Linear=zeros(1,n/2); %vector for interpolated values
V_model_S=zeros(1,n/2); %vector for model values;
V_int_R_Linear=zeros(1,n/2);
V_model_R=zeros(1,n/2);
V_int_I_Linear=zeros(1,n/2);
V_model_I=zeros(1,n/2);

%% Newton's linear interpolation Method %%
for k = 1:2:n
    S_linear(k+1)=S(k)+((S(k+2)-S(k))/((k+2)-k))*((k+1)-k);
    R_linear(k+1)=R(k)+((R(k+2)-R(k))/((k+2)-k))*((k+1)-k);
    I_linear(k+1)=I(k)+((I(k+2)-I(k))/((k+2)-k))*((k+1)-k);
end

for k=1:n/2
    V_model_S(k)=S(2*k);
    V_model_R(k)=R(2*k);
    V_model_I(k)=I(2*k);
    V_int_S_Linear(k)=S_linear(2*k);
    V_int_R_Linear(k)=R_linear(2*k);
    V_int_I_Linear(k)=I_linear(2*k);
end


%%Setting vectors for quadratic interpolation%%
V_int_S_Quadratic=zeros(1,n/2); %vector for interpolated values
V_int_R_Quadratic=zeros(1,n/2);
V_int_I_Quadratic=zeros(1,n/2);

S_Quadratic=S_linear;
R_Quadratic=R_linear;
I_Quadratic=I_linear;

for i=1:(n/2)
    S_Quadratic(2*i)=0;
    R_Quadratic(2*i)=0;
    I_Quadratic(2*i)=0;

end



%% Newton Quadratic Interpolation %%
for k = 2:2:n
    
    if(k < n)
    b0_S = S(k-1);
    b1_S = (S(k+1)-S(k-1))/((k+1)-(k-1));
    b11_S = (S(k+3)-S(k+1))/((k+3)-(k+1));
    b2_S = (b11_S-b1_S)/((k+3)-(k-1));
    S_Quadratic(k)=b0_S+b1_S*(k-(k-1))+b2_S*(k-(k-1))*(k-(k+1));
    b0_R = R(k-1);
    b1_R = (R(k+1)-R(k-1))/((k+1)-(k-1));
    b11_R = (R(k+3)-R(k+1))/((k+3)-(k+1));
    b2_R = (b11_R-b1_R)/((k+3)-(k-1));
    R_Quadratic(k)=b0_R+b1_R*(k-(k-1))+b2_R*(k-(k+1));
    b0_I = I(k-1);
    b1_I = (I(k+1)-I(k-1))/((k+1)-(k-1));
    b11_I = (I(k+3)-I(k+1))/((k+3)-(k+1));
    b2_I = (b11_I-b1_I)/((k+3)-(k-1));
    I_Quadratic(k)=b0_I+b1_I*(k-(k-1))+b2_I*(k-(k+1));

    elseif(k == n)
    b0_S = S(k-3);
    b1_S = (S(k-1)-S(k-3))/((k-1)-(k-3));
    b11_S = (S(k+1)-S(k-1))/((k+1)-(k-1));
    b2_S = (b11_S-b1_S)/((k+3)-(k-1));
    S_Quadratic(k)=b0_S+b1_S*(k-(k-3))+b2_S*(k-(k-3))*(k-(k-1));
    b0_R = R(k-3);
    b1_R = (R(k-1)-R(k-3))/((k-1)-(k-3));
    b11_R = (R(k+1)-R(k-1))/((k+1)-(k-1));
    b2_R = (b11_R-b1_R)/((k+3)-(k-1));
    R_Quadratic(k)=b0_R+b1_R*(k-(k-3))+b2_R*(k-(k-3))*(k-(k-1));
    b0_I = I(k-3);
    b1_I = (I(k-1)-I(k-3))/((k-1)-(k-3));
    b11_I = (I(k+1)-I(k-1))/((k+1)-(k-1));
    b2_I = (b11_I-b1_I)/((k+3)-(k-1));
    I_Quadratic(k)=b0_I+b1_I*(k-(k-3))+b2_I*(k-(k-3))*(k-(k-1));
    
     
    end



   

end

for k=1:n/2
    V_int_S_Quadratic(k)=S_Quadratic(2*k);
    V_int_R_Quadratic(k)=R_Quadratic(2*k);
    V_int_I_Quadratic(k)=I_Quadratic(2*k);
end



%%Error Norm El2 for I(t), R(t), S(t)%%
N_tot = n/2;
Idiff_Linear = V_int_I_Linear-V_model_I;
Rdiff_Linear = V_int_R_Linear-V_model_R;
Sdiff_Linear = V_int_S_Linear-V_model_S;
EL2_I_Linear = sqrt((sum((Idiff_Linear).^2))/(N_tot));
EL2_R_Linear = sqrt((sum((Rdiff_Linear).^2))/(N_tot));
EL2_S_Linear = sqrt((sum((Sdiff_Linear).^2))/(N_tot));
Idiff_Quadratic = V_int_I_Quadratic-V_model_I;
Rdiff_Quadratic = V_int_R_Quadratic-V_model_R;
Sdiff_Quadratic = V_int_S_Quadratic-V_model_S;
EL2_I_Quadratic = sqrt((sum((Idiff_Quadratic).^2))/(N_tot));
EL2_R_Quadratic = sqrt((sum((Rdiff_Quadratic).^2))/(N_tot));
EL2_S_Quadratic = sqrt((sum((Sdiff_Quadratic).^2))/(N_tot));