function output=F1(AA,A_g,h,M1,M2,M3,m1,m2,m3,sigma)

A0=[1 0 0 -M1(1) 0
    0 1 0 -M1(2) 0
    0 0 1 -M2 0
    -1 0 0 m1(1) 0
    0 -1 0 m1(2) 0
    0 0 -1 m2 0
    1 0 0 -m1(1) 0
    0 1 0 -m1(2) 0
    0 0 1 -m2 -1
    -1 0 0 M1(1) 0
    0 -1 0 M1(2) 0
    0 0 -1 M2 1
    0 0 0 (m3+sigma) 0
    0 0 0 (M3-sigma) 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 1
    0 0 0 0 -1];
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
    -AA
    zeros(2,5)];

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
    -AA2
    zeros(2,5)];

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
    -AA3
    zeros(2,5)];
%%
ZZ=zeros(20,5);

output=[A0 ZZ ZZ ZZ
    A_10 A0 ZZ ZZ
    A_20 A_10 A0 ZZ
    A_30 A_20 A_10 A0];
end

