%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2023 by Sichen Yuan
%   Created: 2023/07/07
%   $Revision: 1.0 $  $Date: 2023/07/07 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nodes,nodes_s,theta_e] = ell_cal()

global N_chord
global a_axis
global b_axis

xx0 = ones(N_chord-1,1)*2*pi()/N_chord;
theta = zeros(N_chord,1);
[tempp,fval,exitflag,output] = fsolve(@equal_length_fun,xx0);
for ii = 2:1:N_chord
    theta(ii) = sum(tempp(1:ii-1));
end
theta_e = theta/2/pi();
nodes = [cos(theta) b_axis/a_axis*sin(theta)];

nodes_s = zeros(N_chord,2);
for ii = 1:1:N_chord
    nodes_s(ii,:) = [cos((ii-1)*2*pi()/N_chord) sin((ii-1)*2*pi()/N_chord)];
end
