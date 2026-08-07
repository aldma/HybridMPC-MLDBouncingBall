function output = b(A_g,x,m2,M1,M2,h,sigma,M3,xmax,xmin,m1)


b_0=[zeros(6,1) % // blocks 5-6-7-8
    x           % // block 13
    -m2         % // block 14 why -m2 and not 0 here?
    -x          % // block 15
    M2          % // block 16
    dot(h, x)   % // block 9
    -dot(h, x)  % // block 10
    zeros(4,1)];

b00=[zeros(6,1) % // blocks 5-6-7-8
    -m1     % // block 13
    -m2     % // block 14
    M1   % // block 15
    M2      % // block 16
    sigma   % // block 9
    M3      % // block 10
    zeros(4,1)];

b0=b00+b_0;

%% for one step befor
x=A_g*x;
b_1=[zeros(6,1)
    x
    0
    -x
    0
    dot(h, x)
    -dot(h, x)
    -x+xmax
    x-xmin];

b1=b00+b_1;
%% for two step befor
x=A_g*x;
b_2=[zeros(6,1)
    x
    0
    -x
    0
    dot(h, x)
    -dot(h, x)
    -x+xmax
    x-xmin];

b2=b00+b_2;

%% for two step befor
x=A_g*x;
b_3=[zeros(6,1)
    x
    0
    -x
    0
    dot(h, x)
    -dot(h, x)
    -x+xmax
    x-xmin];

b3=b00+b_3;

output=[b0
    b1
    b2
    b3];
end