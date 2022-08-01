%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yy = equal_length_fun(xx)

global N_chord
global a_axis
global b_axis


nodes = zeros(N_chord,2);
nodes(1,:) = [cos(0) b_axis/a_axis*sin(0)];

for ii=2:1:N_chord
    temm = sum(xx(1:ii-1));
    nodes(ii,:) = [cos(temm) b_axis/a_axis*sin(temm)];
end

yy = zeros(N_chord-1,1);
for ii= 1:1:N_chord-2
    yy(ii) = ((nodes(ii+1,1)-nodes(ii,1))^2+(nodes(ii+1,2)-nodes(ii,2))^2)^0.5-((nodes(ii+2,1)-nodes(ii+1,1))^2+(nodes(ii+2,2)-nodes(ii+1,2))^2)^0.5;
end
yy(N_chord-1) = ((nodes(N_chord,1)-nodes(N_chord-1,1))^2+(nodes(N_chord,2)-nodes(N_chord-1,2))^2)^0.5-((nodes(N_chord,1)-nodes(1,1))^2+(nodes(N_chord,2)-nodes(1,2))^2)^0.5;
