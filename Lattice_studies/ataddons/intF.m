function [F]=intF(x)
u=0.01:0.001:1;
f=(2./u-log(1./u)-2).*exp(-x./u);
F=sum(f)*0.001;

