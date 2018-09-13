function [X] = atdatadp1(RING,d)
% Higher order chromaticities
dp=(-d : d/4 : d);
i=0;
for dd=dp
    i=i+1;
    [~, tunes] = twissring(RING, dd, 1:length(RING)+1);
    nux(i)=tunes(1);nuz(i)=tunes(2);
end

% Get polynome
px = polyfit(dp,nux,3);
pz = polyfit(dp,nuz,3);

% Second order chromaticies
X=[px(3) px(2) px(1)  pz(3) pz(2) pz(1)];

end

