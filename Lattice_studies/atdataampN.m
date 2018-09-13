function [X] = atdataampN(RING)
% Some optics data for atmatch
x=1e-3;

[nux,nuz]=nuampl_alex(RING,[0  x],1,0);
dnuxx=(nux(1,2)-nux(1,1))/x^2; 
dnuzx=(nuz(1,2)-nuz(1,1))/x^2; 

[nux,nuz]=nuampl_alex(RING,[0  x],2,0);
dnuxz=(nux(1,2)-nux(1,1))/x^2; 
dnuzz=(nuz(1,2)-nuz(1,1))/x^2; 

X=[dnuxx dnuzx dnuxz  dnuzz]/1000;

end

