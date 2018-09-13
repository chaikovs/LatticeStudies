

%%
%Lmag = 0.296233 @ 159.2475 A

dL = (0.296233-0.27646)/2;

% SD31          SD 0.7537500E+00 
% SD3S2         SD 0.1250000E+00
% SD0           SD 0.1000000E+00 
% SD3           SD 0.2350000E+00

SD31  = 0.7537500E+00 ;
SD3S2 = 0.1250000E+00;
SD0 = 0.1000000E+00 ;
SD3 =0.2350000E+00;

SD31new = SD31 - dL
SD3S2new = SD3S2 - dL
SD0new = SD0 - dL
SD3new = SD3 - dL

%%

% DI 0.7853980E+00 0.3520000E+00 

Lnom = 0.7853980*0.352;
rho_new = 0.296233/0.7853980

Lnew = 0.7853980*0.37720

Lnew_good = 0.7853980*0.3771756

%%

%Lmag = 0.296351 @ 159.2475 A from tracking

dL_track = (0.296351 - 0.296233)/2;

SD31  =     0.7402085E+00 
SD3S2 =     0.1151135E+00 
SD0     =   0.9011350E-01 
SD3    =    0.2214585E+00 

SD31new = SD31 - dL_track
SD3S2new = SD3S2 - dL_track
SD0new = SD0 - dL_track
SD3new = SD3 - dL_track
