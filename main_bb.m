clear all %#ok
close all
clc

%% flow map parameters
T_s = 0.01;
g = 9.8;
lambda = 0.7;
A_f = [1, T_s; 0, 1];
B_f = [0; 0];
C_f = [-T_s^2*g; -T_s*g];

%% jump map parameters
A_g = [1, -T_s; 0, -lambda];
B_g = [0; 1];
C_g = [0; 0];

%% jump and flow set  parameters
h = [-1, 0];
sigma = 0;

%% feasiblity set
x_1max = 10;
x_1min = -10;
x_2max = 10;
x_2min = -10;
u_max = 0.01;
u_min = -0.01;
xmax=[x_1max; x_2max];
xmin=[x_1min; x_2min];

%% finding max and min of flow map
M11=x_1max;
M12=x_2max;
m11=x_1min;
m12=x_2min;
m1=[m11; m12];
M1=[M11; M12];

%% finding max and min of jump map
M2=u_max;
m2=u_min;

%% max and min on h(x)
cc8 = [h(1);h(2)]>=0;
M3  = h(1)*(cc8(1)*x_1max+(1-cc8(1))*x_1min)+h(2)*(cc8(2)*x_2max+(1-cc8(2))*x_2min);
cc9 = [h(1);h(2)]<=0;
m3  = h(1)*(cc9(1)*x_1max+(1-cc9(1))*x_1min)+h(2)*(cc9(2)*x_2max+(1-cc9(2))*x_2min);

%% calculating AA for x(\ell)=AA*[z(\ell-1),u_f(\ell-1),u(\ell-1)]
AA=AA(A_f,A_g,B_f,B_g,C_f,C_g);

%% calculating A for F1*X<b
A=F1(AA,A_g,h,M11,M12,M2,M3,m11,m12,m2,m3,sigma);

%% cost function parameters
Q_c = 2e-1*eye(2);
Q_d = 2e-1*eye(2);
P   = 1e-1*eye(2);
R_c = 1e-2;
R_d = 1e-2;

%% initial value
x1 = 2;
x2 = 0;
T_f = 349; % simulation time
T_p = 1; % control horizon
N = 4;
optx0 = zeros(5*N,1);

x_1 = nan(1,T_f+1);
x_2 = nan(1,T_f+1);
x_1(1) = x1;
x_2(1) = x2;
u_f = nan(1,T_f);
u_ = nan(1,T_f);

% binary variables
ivar = 4:5:(5*N);

% solver configuration
options = [];

EXT_noptx = 7*N;
EXT_optx0 = zeros(EXT_noptx,1);
EXT_ivar = 4:7:EXT_noptx;
EXT_idxx = sort([6:7:EXT_noptx, 7:7:EXT_noptx]);
EXT_idxa = setdiff(1:EXT_noptx,EXT_idxx);
EXT_A = zeros(size(A,1),EXT_noptx);
EXT_A(:,EXT_idxa) = A;
EXT_Aeq = zeros(2*N,EXT_noptx);
EXT_beq = zeros(2*N,1);
for j = 1:N
    EXT_Aeq((1:2)+2*(j-1),(1:5)+7*(j-1)) = -AA;
    EXT_Aeq((1:2)+2*(j-1),(6:7)+7*(j-1)) = eye(2);
    if j > 1
        EXT_Aeq((1:2)+2*(j-1),(6:7)+7*(j-2)) = -A_g;
        EXT_beq((1:2)+2*(j-1)) = C_g;
    end
end
EXT_S1 = zeros(EXT_noptx,EXT_noptx);
EXT_S2 = zeros(EXT_noptx,1);

for i=1:T_f
    %% b for next steps that x(i) is written as sumation of z_1 and z_2 so we dont need them in b
    x = [x1; x2];

    % optvec = [z1,rho1,u1, z2,rho2,u2, ..., zN,rhoN,uN]
    [S1,S2] = costfunction(AA,Q_c,A_g,Q_d,R_d,R_c,P,x);
    bb = b(A_g,x,x1,x2,m2,M11,M12,M2,h,sigma,M3,xmax,xmin,m11,m12,u_max,u_min);
    
    % solver call
    optsol = miqp(S1,S2,A,bb,[],[],ivar,[],[],optx0,options);

    x_11=AA*optsol(1:5)+A_g*x+C_g;
    x1=x_11(1);
    x2=x_11(2);
    x_1(i+1) = x1;
    x_2(i+1) = x2;
    u_f(i) = optsol(4);
    u_(i) = optsol(5);

    %=======
    % with extended decision vector
    % EXT_optvec = [z1,rho1,u1,x2, z2,rho2,u2,x3, ..., zN,rhoN,uN,xN+1]
    %=======
    EXT_S1(EXT_idxa,EXT_idxa) = S1;
    EXT_S2(EXT_idxa) = S2;
    EXT_b = bb;
    EXT_beq(1:2) = A_g*x + C_g;
    % solver call
    EXT_optsol = miqp(EXT_S1,EXT_S2,EXT_A,EXT_b,EXT_Aeq,EXT_beq,...
        EXT_ivar,[],[],EXT_optx0,options);
    EXT_x_11 = EXT_optsol(6:7);

    % warmstart
    optx0 = optsol;
    EXT_optx0 = EXT_optsol;
end

%% plotting
plotting(T_p,T_f,u_f,x_1,x_2,u_)