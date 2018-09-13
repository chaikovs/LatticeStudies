function [ Z] = matrix_coin(s, X0 )
%MATRIX_COIN Summary of this function goes here
%   Detailed explanation goes here
l=length(s);
ds=(s(l)-s(1))/l;
X0(4:6)=X0(4:6)/2.99792458e8;
beta=0;
rho=5.36;
gap=0.038;
K=0.17;
phi=K*gap/rho*(1+sin(beta))^2/cos(beta);

D=[1,ds;0 1];
C=[1,0;-tan(beta-phi)/rho 1];
Z(:,1)=[X0(2)  X0(5)];
coin=0;
for i=2:l
    Z(:,i)=D*Z(:,i-1);
    if (s(i)>0.505) && (coin==0)
        coin=1;
        Z(:,i)=C*Z(:,i-1);
    end
end

