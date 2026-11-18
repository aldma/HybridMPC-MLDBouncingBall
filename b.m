function output=b(A_g,x,x1,x2,m2,M11,M12,M2,h,sigma,M3,xmax,xmin,m11,m12,u_max,u_min)


 b_0=[zeros(6,1)
    x1
    x2
    0
    -x1
    -x2
    0
    0+h(1,1)*x1+h(1,2)*x2
    -h(1,1)*x1-h(1,2)*x2
    0
    0
    0
    0
    0
    0];

 b00=[zeros(6,1)
    -m11
    -m12
    -m2
    +M11
    +M12
    M2
    0+sigma
    M3
    0
    0
    0
    0
    u_max
    -u_min];
b0=b00+b_0;

%% for one step befor
x=A_g*x;
b_1=[zeros(6,1)
    x
    0
    -x
    0
    h(1,1)*x(1)+h(1,2)*x(2)
    -h(1,1)*x(1)-h(1,2)*x(2)
    -x+xmax
    x-xmin
    u_max
    -u_min];

b1=b00+b_1;
%% for two step befor
x=A_g*x;
b_2=[zeros(6,1)
    x
    0
    -x
    0
    h(1,1)*x(1)+h(1,2)*x(2)
    -h(1,1)*x(1)-h(1,2)*x(2)
    -x+xmax
    x-xmin
    u_max
    -u_min];

b2=b00+b_2;

%% for two step befor
x=A_g*x;
b_3=[zeros(6,1)
    x
    0
    -x
    0
    h(1,1)*x(1)+h(1,2)*x(2)
    -h(1,1)*x(1)-h(1,2)*x(2)
    -x+xmax
    x-xmin
    u_max
    -u_min];

b3=b00+b_3;
output=[b0
    b1
    b2
    b3];
end