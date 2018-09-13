function [tx, ty] = tracem44_AT2(ring, dp)
%  [tx, ty] = tracem44(dp)

tx = zeros(1,length(dp));
ty = zeros(1,length(dp));

for i =1:length(dp)   
    M44 = findm44(ring,dp(i));
    MX = M44([1 2 ], [1 2]);
    MY = M44([3 4], [ 3 4]);
    tx(i) = trace(MX);
    ty(i) = trace(MY);
end