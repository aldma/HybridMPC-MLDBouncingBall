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
Adyn = AA(A_f,A_g,B_f,B_g,C_f,C_g);

%% calculating A for F1*X<b
A = F1(Adyn,A_g,h,M1,M2,M3,m1,m2,m3,sigma);

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

% lower and upper bounds
lb = -inf(5*N,1);
ub = inf(5*N,1);
lb(5:5:end) = u_min;
ub(5:5:end) = u_max;

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
[EXT_Aeq, EXT_beq] = make_Aeq_beq(N,Adyn,A_g,C_g);
EXT_S1 = zeros(EXT_noptx,EXT_noptx);
EXT_S2 = zeros(EXT_noptx,1);
EXT_lb = -inf(EXT_noptx,1);
EXT_ub = inf(EXT_noptx,1);
EXT_lb(EXT_idxa) = lb;
EXT_ub(EXT_idxa) = ub;
EXT_lb(EXT_idxx) = repmat(xmin,N,1);
EXT_ub(EXT_idxx) = repmat(xmax,N,1);

EXT_x1 = x1;
EXT_x2 = x2;
EXT_x_1 = nan(1,T_f+1);
EXT_x_2 = nan(1,T_f+1);
EXT_x_1(1) = EXT_x1;
EXT_x_2(1) = EXT_x2;
EXT_u_f = nan(1,T_f);
EXT_u_ = nan(1,T_f);

for i = 1:T_f
    % current state
    x = [x1; x2];
    EXT_x = [EXT_x1; EXT_x2];

    %============================================================%
    % compact decision vector: implicit state x
    % optvec = [z1,rho1,u1, z2,rho2,u2, ..., zN,rhoN,uN]
    %============================================================%
    [S1,S2] = costfunction(Adyn,Q_c,A_g,Q_d,R_d,R_c,P,x);
    bb = b(A_g,x,m2,M1,M2,h,sigma,M3,xmax,xmin,m1);
    % solver call
    optsol = miqp(S1,S2,A,bb,[],[],ivar,lb,ub,optx0,options);
    % store values
    x_11 = Adyn*optsol(1:5)+A_g*x+C_g;
    x_1(i+1) = x_11(1);
    x_2(i+1) = x_11(2);
    u_f(i) = optsol(4);
    u_(i) = optsol(5);

    
    %============================================================%
    % extended decision vector: explicit state x
    % EXT_optvec = [z1,rho1,u1,x2, z2,rho2,u2,x3, ..., zN,rhoN,uN,xN+1]
    %============================================================%
    [tmpS1,tmpS2] = costfunction(Adyn,Q_c,A_g,Q_d,R_d,R_c,P,EXT_x);
    EXT_S1(EXT_idxa,EXT_idxa) = tmpS1;
    EXT_S2(EXT_idxa) = tmpS2;
    EXT_beq = update_beq(EXT_beq,A_g,C_g,EXT_x);
    % solver call
    EXT_optsol = miqp(EXT_S1,EXT_S2,EXT_A,EXT_b,EXT_Aeq,EXT_beq,...
        EXT_ivar,EXT_lb,EXT_ub,EXT_optx0,options);
    % store values
    EXT_u_f(i) = EXT_optsol(4);
    EXT_u_(i) = EXT_optsol(5);
    EXT_x_1(i+1) = EXT_optsol(6);
    EXT_x_2(i+1) = EXT_optsol(7);

    % warmstart
    optx0 = optsol;
    EXT_optx0 = EXT_optsol;

    % next state
    x1 = x_11(1);
    x2 = x_11(2);
    EXT_x1 = EXT_optsol(6);
    EXT_x2 = EXT_optsol(7);
end

%% plotting
figure
plotting(T_p,T_f,u_f,x_1,x_2,u_)

figure
plotting(T_p,T_f,EXT_u_f,EXT_x_1,EXT_x_2,EXT_u_)
drawnow

%% auxiliary functions for extended MIQP formulation
function [Aeq, beq] = make_Aeq_beq(N,Adyn,A_g,C_g)
%MAKE_AEQ_BEQ
% see eq (23), paper 280.pdf
    noptx = 7*N;
    Aeq = zeros(2*N,noptx);
    beq = zeros(2*N,1);
    for j = 1:N
        Aeq((1:2)+2*(j-1),(1:5)+7*(j-1)) = -Adyn;
        Aeq((1:2)+2*(j-1),(6:7)+7*(j-1)) = eye(2);
        if j > 1
            Aeq((1:2)+2*(j-1),(6:7)+7*(j-2)) = -A_g;
            beq((1:2)+2*(j-1)) = C_g;
        end
    end
end

function beq = update_beq(beq,A_g,C_g,xnow)
%UPDATE_BEQ
% see eq (23), paper 280.pdf
    beq(1:2) = A_g*xnow + C_g;
end