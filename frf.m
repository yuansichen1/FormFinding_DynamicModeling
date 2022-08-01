%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,G]=frf(X)
[U,S,V] = svd(X);
j=1;
for i=1:1:length(diag(S))
    if abs(S(i,i))>1e-8
        indx(j)=i;
        j=j+1;
    end
end
indx = indx';
F = U(:,indx);
G = S(indx,indx)*V(:,indx)';      