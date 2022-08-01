%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function setup_struc_parameter(End_Spec, Node_Spec, BC_Spec)

global Num_Member           % number of members
global Num_Node             % number of nodes
global Member_Ends          % node numbers (Ii, Ji) of members
global Member_Length_des    % length (Li) of members
global Node_Desired         % global coordinates (Xi, Yi, Zi) of nodes
global BCs                  % info on boundar conditions
global Num_BCNode           % number of boundary nodes

% Check number of members and number of nodes
Num_Member = length(End_Spec(:,1));
if Num_Member == 1
    error('Number of members = 1. This Toolbox is not valid for single-member truss.')
end
Num_Node = length(Node_Spec(:,1));
Member_Ends = End_Spec;
Node_Desired = Node_Spec;
BCs = BC_Spec;
[Num_BCNode,ttteeee] = size(BC_Spec);

if Num_Node > 200
    error('The Toolbox is limited up to 200 nodes in one reflector surface.')
end

% Compute the length, direction of members
Member_Length_des = zeros(Num_Member,1);
Member_DirCosin = zeros(Num_Member, 3);
for i = 1: Num_Member
    i0 = Member_Ends(i,1);
    iL = Member_Ends(i, 2);
    dx = Node_Desired(iL, 1) - Node_Desired(i0, 1);
    dy = Node_Desired(iL, 2) - Node_Desired(i0, 2);
    dz = Node_Desired(iL, 3) - Node_Desired(i0, 3);
    Member_Length_des(i) = sqrt(dx^2 + dy^2 + dz^2);
    if Member_Length_des(i) == 0, 
        error('Warning: the length of a member = 0!'); 
    end
    Member_DirCosin(i,:) = [dx dy dz]/Member_Length_des(i);
end