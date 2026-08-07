function output=F1(AA,A_g,h,M1,M2,M3,m1,m2,m3,sigma)

A0=[1 0 0 -M1(1) 0  % z1 <= M1*rho // block 5, eq (22), paper 280
    0 1 0 -M1(2) 0
    0 0 1 -M2 0     % z2 <= M2*rho // block 6
    -1 0 0 m1(1) 0  % z1 >= m1*rho // block 7
    0 -1 0 m1(2) 0
    0 0 -1 m2 0     % z2 >= m2*rho // block 8
    1 0 0 -m1(1) 0  % z1 - xhat + m1*(1-rho) <= 0 // block 13
    0 1 0 -m1(2) 0
    0 0 1 -m2 -1    % z2 - u + m2*(1-rho) <= 0 // block 14 
    -1 0 0 M1(1) 0  % xhat - z1 - M1*(1-rho) <= 0 // block 15
    0 -1 0 M1(2) 0
    0 0 -1 M2 1     % u - z2 - M2*(1-rho) <= 0 // block 16
    0 0 0 (m3+sigma) 0 % // block 9
    0 0 0 (M3-sigma) 0 % // block 10
    zeros(4,5)];

%% for one step befor
H1=h(1,1)*AA(1,:)+h(1,2)*AA(2,:);

A_10=[zeros(6,5)
    -AA
    zeros(1,5)
    AA
    zeros(1,5)
    -H1
    H1
    AA
    -AA];

%% for two step befor
AA2=A_g*AA;

H1=h(1,1)*AA2(1,:)+h(1,2)*AA2(2,:);

A_20=[zeros(6,5)
    -AA2
    zeros(1,5)
    AA2
    zeros(1,5)
    -H1
    H1
    AA2
    -AA2];

%% for three step befor
AA3=A_g^2*AA;

H1=h(1,1)*AA3(1,:)+h(1,2)*AA3(2,:);

A_30=[zeros(6,5)
    -AA3
    zeros(1,5)
    AA3
    zeros(1,5)
    -H1
    H1
    AA3
    -AA3];

%%
ZZ=zeros(18,5);

output=[A0 ZZ ZZ ZZ
    A_10 A0 ZZ ZZ
    A_20 A_10 A0 ZZ
    A_30 A_20 A_10 A0];
end

