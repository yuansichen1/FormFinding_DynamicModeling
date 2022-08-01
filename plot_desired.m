%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_desired(display_node_no)

%	Plot configurations of truss
%       display_node_no == 0    do NOT show node and memeber numbers
%       display_node_no == 1    show node numbers
%       display_node_no == 2    show memeber numbers
%       display_node_no == 3    show node and memeber numbers

global Num_Member           % number of members
global Num_Node             % number of nodes
global Member_Ends          % node numbers (Ii, Ji) of members
global Member_Length_des    % length (Li) of members
global Node_Desired         % global coordinates (Xi, Yi, Zi) of nodes
global BCs                  % info on boundar conditions
global Num_BCNode           % number of boundary nodes

if Num_Node > 200
    error('The Toolbox is limited up to 200 nodes in one reflector surface.')
end

if nargin < 1 || isempty(display_node_no), display_node_no = 1; end   % default: do not show node and memeber numbers
clf

% Figure boundaries
xL_node = max(Node_Desired(:,1)); x0_node = min(Node_Desired(:,1));
yL_node = max(Node_Desired(:,2)); y0_node = min(Node_Desired(:,2));
zL_node = max(Node_Desired(:,3)); z0_node = min(Node_Desired(:,3));
ee = 0.1;
xmax = xL_node + ee*(xL_node-x0_node);
xmin = x0_node - ee*(xL_node-x0_node);
ymax = yL_node + ee*(yL_node-y0_node);
ymin = y0_node - ee*(yL_node-y0_node);
zmax = zL_node + ee*(zL_node-z0_node);
zmin = z0_node - ee*(zL_node-z0_node);

% Parametera for locating node numbers
x_shift = max(abs(Node_Desired(:,1)))/30;
y_shift = max(abs(Node_Desired(:,2)))/30;
z_shift = max(abs(Node_Desired(:,3)))/30;

disp(' See figure for the configuration of the truss')
disp(' where star (*)  for a ball-and-socket joint')
disp('       triangle (up) for a slotted roller joint;')
disp('       circle for plane roller joint or a short link')
disp(' ')

% Start plotting members
j0 = Member_Ends(1, 1); jL = Member_Ends(1, 2);
x = [Node_Desired(j0,1) Node_Desired(jL,1)];
y = [Node_Desired(j0,2) Node_Desired(jL,2)];
z = [Node_Desired(j0,3) Node_Desired(jL,3)];
plot3(x,y,z)
if (display_node_no == 2) || (display_node_no == 3) % add member number
    tx_shft = Member_Length_des(1)/20;
    text((x(2)+x(1))/2+tx_shft, (y(2)+y(1))/2-tx_shft, (z(2)+z(1))/2, '(1)');  
end
axis([xmin xmax ymin ymax zmin zmax])
xlabel('x'); ylabel('y'); zlabel('z')
hold
for ib = 2:Num_Member
    j0 = Member_Ends(ib, 1); jL = Member_Ends(ib, 2);
    x = [Node_Desired(j0,1) Node_Desired(jL,1)];
    y = [Node_Desired(j0,2) Node_Desired(jL,2)];
    z = [Node_Desired(j0,3) Node_Desired(jL,3)];
    plot3(x,y, z)
    if (display_node_no == 2) || (display_node_no == 3)   % show member numbers
        tx_shft = Member_Length_des(ib)/20;
        text((x(2)+x(1))/2+tx_shft, (y(2)+y(1))/2-tx_shft, (z(2)+z(1))/2, ['(',num2str(ib),')']);  % add member number
    end
end
    
% Start plotting nodes
node_indx = ones(1, Num_Node);
% mark boundary nodes
for i = 1:Num_BCNode
    node = BCs(i,1);
    node_indx(node) = 0;
    x = Node_Desired(node,1);
    y = Node_Desired(node,2);
    z = Node_Desired(node,3);
    BC_Type = BCs(i,2);
    if BC_Type == 1
        plot3(x,y,z, 'r*')
    elseif BC_Type == 2 
        plot3(x,y,z, 'rx')
    elseif BC_Type == 3
        plot3(x,y,z, 'r^')
    else
        plot3(x,y,z, 'ro')
    end
    if (display_node_no == 1) || (display_node_no == 3)  % show node numbers
            text(x+x_shift, y-y_shift, z+z_shift, num2str(node));  % add node number
    end
end

% plot rest of the nodes
for i = 1:Num_Node
    if node_indx(i) == 1
        x = Node_Desired(i,1);
        y = Node_Desired(i,2);
        z = Node_Desired(i,3);
        plot3(x,y, z,'k.')
        if (display_node_no == 1) || (display_node_no == 3)   % show node number
            text(x+x_shift, y-y_shift, z-z_shift, num2str(i));  % add node number
        end
    end
end
% title('Desired configuration')

hold
disp(' ')