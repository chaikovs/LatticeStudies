
tx=thomx;
%%

%The following will work fine
findorbit6(tx);
ringpass(tx(:),[0;0;0;0;0;0],10,'reuse')

%%
%The following produces a crash in Matlab R2016a
findorbit6(tx(:));
ringpass(tx,[0;0;0;0;0;0],10,'reuse')