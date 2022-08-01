%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F_p, D_p, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b, ratio_zeta, w_b)

% n_r is the number of rings
% n_s is the number of subdivision
% n_c is the converging ring number
% rou is the power rouficient rou
% c_b is the connections of a node at boundary layer
% c_b = 0 when no boundary treatment
% w_b = 0 when entire surface is the working surface

global F_ref
global D_ref
global F_ref_eff
global D_ref_eff
global D_circle
global D_circle_eff
global F_parent
global D_parent
global e_offset
global e_eff
global F_p_eff
global D_p_eff
global n_r_bar
global flag_effect
global flag_shape
global flag_update

if flag_update == 0
    if flag_effect == 0
        ratio_off = e_offset/D_parent*2;
    else
        n_i_i_temp = 0;
        if n_c == n_r
            n_i_i_temp = n_r-1;
        else        
            for i_j_k = 1:1:n_r-1
                if i_j_k < n_c
                    n_i_i_temp_temp = (1+(ratio_zeta-1)/2);
                else
                    n_i_i_temp_temp = (1+((i_j_k-n_c)/(n_r-n_c-1))^(rou)*(1-ratio_zeta)+(ratio_zeta-1)/2);
                end
                n_i_i_temp = n_i_i_temp_temp + n_i_i_temp;
            end
        end

        D_ref = 4*F_ref_eff*sin((n_i_i_temp+w_b)/n_i_i_temp*asin(D_ref_eff/4/F_ref_eff));
        F_ref = F_p_eff;
        
        ratio_off_eff = e_eff/D_p_eff*2;
        angle_eff = atan((1+ratio_off_eff)*D_p_eff/8/F_p_eff);
        angle_off = angle_eff;

        D_parent = D_p_eff+(D_ref-D_ref_eff)*cos(angle_off);
        F_parent = F_p_eff;
        e_offset = e_eff-(D_ref-D_ref_eff)*cos(angle_off)/2;

        ratio_off = e_offset/D_parent*2;
    end
else
    ratio_off = e_offset/D_parent*2;
end
D_p = D_parent;
F_p = F_parent;
e_p = e_offset;


%% Step 1 Calculate the geometry and Set up the references
e_offset = ratio_off*D_parent/2;
R_c = 1/2*(D_parent/2-e_offset);

H_parent = D_parent^2/16/F_parent;
e_zz = ratio_off^2*D_parent^2/16/F_parent;
z_o = -1/2*e_zz-1/2*H_parent;
x_o = e_offset+R_c;
% y_o = 0;

angle_off = atan((1+ratio_off)*D_parent/8/F_parent);
a_ellipse = R_c/cos(angle_off);
b_ellipse = R_c;

bb2 = sin(angle_off)^2/4/F_parent;
bb1 = sin(angle_off)*x_o/2/F_parent+cos(angle_off);
bb0 = z_o+x_o^2/4/F_parent;
H_p = (-bb1+(bb1^2-4*bb2*bb0)^0.5)/2/bb2;

F = F_parent;
D = 2*a_ellipse;


R_s = 2*F;
H = R_s-(R_s^2-D^2/4)^0.5;

a_ellipsoid = R_s;
b_ellipsoid = 4*F/D*R_c;
c_ellipsoid = R_s;
z_prime = H-R_s;
 
theta = acos((R_s-H)/R_s);

n_i_i_temp = 0;
n_i_i = zeros(n_r,1);
if flag_effect == 0
    if n_c == n_r
        n_i_i = [1:1:n_r]';
    else
        for i_j_k = 1:1:n_r
            if i_j_k < n_c
                n_i_i_temp_temp = (1+(ratio_zeta-1)/2);
            else
                n_i_i_temp_temp = (1+((i_j_k-n_c)/(n_r-n_c))^(rou)*(1-ratio_zeta)+(ratio_zeta-1)/2);
            end
            n_i_i(i_j_k) = n_i_i_temp_temp + n_i_i_temp;
            n_i_i_temp = n_i_i(i_j_k);
        end
    end
else
    if n_c == n_r
        n_i_i = [1:1:n_r]';
        n_i_i(n_r) = n_r-1+w_b;
    else        
        for i_j_k = 1:1:n_r-1
            if i_j_k < n_c
                n_i_i_temp_temp = (1+(ratio_zeta-1)/2);
            else
                n_i_i_temp_temp = (1+((i_j_k-n_c)/(n_r-n_c-1))^(rou)*(1-ratio_zeta)+(ratio_zeta-1)/2);
            end
            n_i_i(i_j_k) = n_i_i_temp_temp + n_i_i_temp;
            n_i_i_temp = n_i_i(i_j_k);
        end
        n_i_i(n_r) = n_i_i(n_r-1) + w_b;
    end
end
theta_i = zeros(n_r,1);
for i_j_k =1:1:n_r
    theta_i(i_j_k) = theta/n_i_i(n_r)*n_i_i(i_j_k);
end
%% Step 2 Nodes before convergence
if c_b == 0
    n_t = (1+n_r)/2*n_r*n_s+1;
else
    n_t = round((n_r/2+2/(c_b+1))*(n_r-1)*n_s+1);
end
if c_b == 0
    n_i_rim = n_s*n_r;
else
    n_i_rim = 2*(n_r-1)*n_s/(c_b+1);
end
if c_b == 0
    m_t = (3*n_r+1)/2*n_r*n_s;
else
    m_t = (3*n_r+2)/2*(n_r-1)*n_s;
end
Indx_node_load = [1:1:(n_t-n_i_rim)]';
if flag_effect == 0
    m_nobc = m_t;
else
    if c_b == 0
        m_nobc = m_t-(3*n_r-1)*n_s;
    else
        m_nobc = m_t-2*(n_r-1)*n_s;
    end
end

B_C = [[(n_t-n_i_rim+1):1:n_t]',ones(n_i_rim,1),zeros(n_i_rim,3)];

Node = zeros(n_t,5);
Node(1,:) = [0,1,0,0,H];
k = 2;
for i = 1:1:(n_c-1)
    n_i = i*n_s;
    xij = R_s*sin(theta_i(i));
    yij = 0;
    zij = R_s*cos(theta_i(i))-R_s+H;
    xiji = R_s*sin(theta_i(i))*cos(2*pi()/n_s);
    yiji = R_s*sin(theta_i(i))*sin(2*pi()/n_s);
    NNi = ((xij-xiji)^2+(yij-yiji)^2)^0.5;
    oNp = (xij^2+yij^2)^0.5;
    gamma = acos(NNi/2/R_s);
    oa = R_s*sin(gamma);
    beta = pi()-2*gamma;
    theta_p = zeros(i-1,1);
    phi_p = zeros(i-1,1);
    for p=1:1:(i-1)
        beta_p = p*beta/i;
        ob_p = oa/sin(gamma+beta_p);
        bpNp = NNi/2-ob_p*cos(gamma+beta_p);
        bbp = zij+R_s-H;
        theta_p(p) = acos(bbp/ob_p);
        obp = ob_p*sin(theta_p(p));
        phi_p(p) = acos((obp^2+oNp^2-bpNp^2)/(2*obp*oNp));
    end
    for q = 1:1:n_s
        j_bar = 1+(q-1)*i;
        xij_bar = R_s*sin(theta_i(i))*cos((q-1)*2*pi()/n_s);
        yij_bar = R_s*sin(theta_i(i))*sin((q-1)*2*pi()/n_s);
        zij_bar = R_s*cos(theta_i(i))-R_s+H;
        Node(k,:) = [i,j_bar,xij_bar,yij_bar,zij_bar];
        k = k+1;
        if i>1
            alpha_p = phi_p+(q-1)*2*pi()/n_s;
            for p = 1:1:(i-1)
                xij_bar = R_s*sin(theta_p(p))*cos(alpha_p(p));
                yij_bar = R_s*sin(theta_p(p))*sin(alpha_p(p));
                zij_bar = R_s*cos(theta_p(p))-R_s+H;
                Node(k,:) = [i,(j_bar+p),xij_bar,yij_bar,zij_bar];
                k = k+1;
            end
        end
    end
end
%% Step 3 Nodes within convergence
if n_c < n_r
    n_i = n_c*n_s;
    xij = R_s*sin(theta_i(n_c));
    yij = 0;
    zij = R_s*cos(theta_i(n_c))-R_s+H;
    xiji = R_s*sin(theta_i(n_c))*cos(2*pi()/n_s);
    yiji = R_s*sin(theta_i(n_c))*sin(2*pi()/n_s);
    ApNp = (xij^2+yij^2)^0.5;
    NNi = ((xij-xiji)^2+(yij-yiji)^2)^0.5;
    oNp = (xij^2+yij^2)^0.5;
    gamma = acos(NNi/2/R_s);
    oa = R_s*sin(gamma);
    beta = pi()-2*gamma;
    beta_p = beta/2;
    ob_p = oa/sin(gamma+beta_p);
    bpNp = NNi/2-ob_p*cos(gamma+beta_p);
    bbp = zij+R_s-H;
    theta_pp = acos(bbp/ob_p);
    obp = ob_p*sin(theta_pp);
    alpha_p = acos((obp^2+oNp^2-bpNp^2)/(2*obp*oNp));
    xnc = R_s*sin(theta_pp)*cos(alpha_p);
    ync = R_s*sin(theta_pp)*sin(alpha_p);
    zeta1 = (yiji-yij)/(xiji-xij);
    zeta2 = yij-zeta1*xij;
    Hnc = abs(zeta1*xnc-ync+zeta2)/(zeta1^2+1)^0.5;
    r_nc = ((NNi/2)^2+Hnc^2)/2/Hnc;

    for i = n_c:1:(n_r-1)
        n_i = i*n_s;
        xij = R_s*sin(theta_i(i));
        yij = 0;
        xiji = R_s*sin(theta_i(i))*cos(2*pi()/n_s);
        yiji = R_s*sin(theta_i(i))*sin(2*pi()/n_s);
        ApNp = (xij^2+yij^2)^0.5;

        if n_c == n_r - 1
            r_i = r_nc;
        else
            if flag_effect == 0
                r_i = 1/((1/(R_s*sin(theta_i(i)))-1/(r_nc))*((i-n_c)/(n_r-n_c))^rou+1/r_nc);
            else
                r_i = 1/((1/(R_s*sin(theta_i(i)))-1/(r_nc))*((i-n_c)/(n_r-1-n_c))^rou+1/r_nc);
            end
        end

        NpNpi = ((xij-xiji)^2+(yij-yiji)^2)^0.5;
        apNp = NpNpi/2;    
        mu = asin(apNp/r_i)*2;
        phi = asin(apNp/ApNp);
        ApNp_p = zeros(i-1,1);
        nu_p = zeros(i-1,1);
        for p=1:1:(i-1)
            mu_p = p*mu/i;
            NpNp_p = (2*r_i^2-2*r_i^2*cos(mu_p))^0.5;
            tao = phi-mu/2;
            ApNp_p(p) = (NpNp_p^2+ApNp^2-2*NpNp_p*ApNp*cos((pi()-mu_p)/2-tao))^0.5;
            nu_p(p) = acos((ApNp_p(p)^2+ApNp^2-NpNp_p^2)/2/ApNp_p(p)/ApNp);
        end
        for q = 1:1:n_s
            j_bar = 1+(q-1)*i;
            xij_bar = R_s*sin(theta_i(i))*cos((q-1)*2*pi()/n_s);
            yij_bar = R_s*sin(theta_i(i))*sin((q-1)*2*pi()/n_s);
            zij_bar = R_s*cos(theta_i(i))-R_s+H;
            Node(k,:) = [i,j_bar,xij_bar,yij_bar,zij_bar];
            k = k+1;
            if i>1
                alpha_p = nu_p+(q-1)*2*pi()/n_s;
                for p = 1:1:(i-1)
                    xij_bar = ApNp_p(p)*cos(alpha_p(p));
                    yij_bar = ApNp_p(p)*sin(alpha_p(p));
                    zij_bar = (R_s^2-xij_bar^2-yij_bar^2)^0.5-R_s+H;
                    Node(k,:) = [i,(j_bar+p),xij_bar,yij_bar,zij_bar];
                    k = k+1;
                end
            end
        end
    end
end
%% Step 4 Nodes on the aperture
if c_b == 0
    for j = 1:1:n_r*n_s
        xij_bar = D/2*cos((j-1)*2*pi()/n_r/n_s);
        yij_bar = D/2*sin((j-1)*2*pi()/n_r/n_s);
        zij_bar = 0;
        Node(k,:) = [n_r,j,xij_bar,yij_bar,zij_bar];
        k = k+1;
    end
else
    for j = 1:1:(n_r-1)*n_s/((c_b+1)/2)
        xij_bar = D/2*cos((j-1)*2*pi()/((n_r-1)*n_s/((c_b+1)/2)));
        yij_bar = D/2*sin((j-1)*2*pi()/((n_r-1)*n_s/((c_b+1)/2)));
        zij_bar = 0;
        Node(k,:) = [n_r,j,xij_bar,yij_bar,zij_bar];
        k = k+1;
    end
end       
%% Step 5 Member indices
% Member is counted by layers
M_indx = zeros(m_t,2);
index1 = findnodeindex(Node(:,1:2),[0,1]);
for jj=1:1:n_s
    index2 = findnodeindex(Node(:,1:2),[1,jj]);
    M_indx(jj,:) = [index1,index2];
end
k = n_s+1;
n_i = 1*n_s;
for j = 1:1:(n_i-1)
    index1 = findnodeindex(Node(:,1:2),[1,j]);
    index2 = findnodeindex(Node(:,1:2),[1,j+1]);
    M_indx(k,:) = [index1,index2];
    k = k+1;
end
index1 = findnodeindex(Node(:,1:2),[1,n_i]);
index2 = findnodeindex(Node(:,1:2),[1,1]);
M_indx(k,:) = [index1,index2];
k = k+1;

if c_b == 0
    for i = 2:1:n_r
        n_i = (i-1)*n_s;
        n_ip = i*n_s;
        index1 = findnodeindex(Node(:,1:2),[i-1,1]);
        index2 = findnodeindex(Node(:,1:2),[i,n_ip]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
        index2 = findnodeindex(Node(:,1:2),[i,1]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
        index2 = findnodeindex(Node(:,1:2),[i,2]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
        i_start = 2;
        for j = 2:1:n_i
            in_onsub = ([1:1:n_s]-1).*(i-1)+1;
            if min(abs(j - in_onsub)) < 1e-6
                index1 = findnodeindex(Node(:,1:2),[i-1,j]);
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
                i_start = i_start +1;
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
                i_start = i_start +1;
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
            else
                index1 = findnodeindex(Node(:,1:2),[i-1,j]);
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
                i_start = i_start +1;
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
            end
        end

        n_i = i*n_s;
        for j = 1:1:(n_i-1)
            index1 = findnodeindex(Node(:,1:2),[i,j]);
            index2 = findnodeindex(Node(:,1:2),[i,j+1]);
            M_indx(k,:) = [index1,index2];
            k = k+1;
        end
        index1 = findnodeindex(Node(:,1:2),[i,n_i]);
        index2 = findnodeindex(Node(:,1:2),[i,1]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
    end
else
    for i = 2:1:(n_r-1)
        n_i = (i-1)*n_s;
        n_ip = i*n_s;
        index1 = findnodeindex(Node(:,1:2),[i-1,1]);
        index2 = findnodeindex(Node(:,1:2),[i,n_ip]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
        index2 = findnodeindex(Node(:,1:2),[i,1]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
        index2 = findnodeindex(Node(:,1:2),[i,2]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
        i_start = 2;
        for j = 2:1:n_i
            in_onsub = ([1:1:n_s]-1).*(i-1)+1;
            if min(abs(j - in_onsub)) < 1e-6
                index1 = findnodeindex(Node(:,1:2),[i-1,j]);
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
                i_start = i_start +1;
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
                i_start = i_start +1;
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
            else
                index1 = findnodeindex(Node(:,1:2),[i-1,j]);
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
                i_start = i_start +1;
                index2 = findnodeindex(Node(:,1:2),[i,i_start]);
                M_indx(k,:) = [index1,index2];
                k = k+1;
            end
        end

        n_i = i*n_s;
        for j = 1:1:(n_i-1)
            index1 = findnodeindex(Node(:,1:2),[i,j]);
            index2 = findnodeindex(Node(:,1:2),[i,j+1]);
            M_indx(k,:) = [index1,index2];
            k = k+1;
        end
        index1 = findnodeindex(Node(:,1:2),[i,n_i]);
        index2 = findnodeindex(Node(:,1:2),[i,1]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
    end
    
    index1 = findnodeindex(Node(:,1:2),[n_r,1]);
    index2 = findnodeindex(Node(:,1:2),[n_r-1,1]);
    M_indx(k,:) = [index1,index2];
    k = k+1; 
    for ijij = 1:1:(c_b-1)/2
        index2 = findnodeindex(Node(:,1:2),[n_r-1,((n_r-1)*n_s-ijij+1)]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
        index2 = findnodeindex(Node(:,1:2),[n_r-1,ijij+1]);
        M_indx(k,:) = [index1,index2];
        k = k+1;        
    end
  
    n_i = 2*(n_r-1)*n_s/(c_b+1);
    for j = 2:1:n_i
        for ijij = 1:1:c_b
            index1 = findnodeindex(Node(:,1:2),[n_r,j]);
            index2 = findnodeindex(Node(:,1:2),[n_r-1,(c_b+1)/2*(j-1)-(c_b-1)/2+ijij]);
            M_indx(k,:) = [index1,index2];
            k = k+1;   
        end
    end
    for j = 1:1:(n_i-1)
        index1 = findnodeindex(Node(:,1:2),[n_r,j]);
        index2 = findnodeindex(Node(:,1:2),[n_r,j+1]);
        M_indx(k,:) = [index1,index2];
        k = k+1;
    end
    index1 = findnodeindex(Node(:,1:2),[n_r,n_i]);
    index2 = findnodeindex(Node(:,1:2),[n_r,1]);
    M_indx(k,:) = [index1,index2];
    k = k+1;
end    
%% mapping
global N_chord
global a_axis
global b_axis

a_axis = a_ellipse;
b_axis = b_ellipse;

Node_e_x = zeros(n_t,1);
Node_e_y = zeros(n_t,1);
Node_e_z = ones(n_t,1)*H;
a_temp = ((1-(Node(:,5)-z_prime).^2/a_ellipsoid^2)*a_ellipsoid^2).^0.5;


if n_s == 4
    k_cal_cr = 10;
    i_cal_cr = 13;
else if n_s == 6
        k_cal_cr = 7;
        i_cal_cr = 9;
    else
        k_cal_cr = 5;
        i_cal_cr = 7;        
    end
end

for i_nc_now = 1:1:n_r-1
    nodal_indx_start_onring = (1+i_nc_now-1)/2*(i_nc_now-1)*n_s+2;
    nodal_indx_end_onring = (1+i_nc_now)/2*i_nc_now*n_s+1;
    
    if i_nc_now < k_cal_cr+1
        N_chord = (i_nc_now)*n_s;
        [nodes_e,nodes_s,theta_e] = ell_cal();
        indp_va_rim = ([1:1:N_chord]'-1)/N_chord;
        fun_val_rim = theta_e;
        [fitobject_rim,gof_rim] = fit([indp_va_rim;1],[fun_val_rim;1],'cubicspline');
        theta_rotation = atan2(Node(nodal_indx_start_onring:nodal_indx_end_onring,4),Node(nodal_indx_start_onring:nodal_indx_end_onring,3));
        for kkjj = 1:1:length(theta_rotation)
            if theta_rotation(kkjj)+2*pi() < 2*pi()
                theta_rotation(kkjj) = theta_rotation(kkjj)+2*pi();
            end
        end
        theta_ellipse = feval(fitobject_rim,theta_rotation/2/pi())*2*pi();
        Node_e_x(nodal_indx_start_onring:nodal_indx_end_onring,1) = cos(theta_ellipse).*a_temp(nodal_indx_start_onring:nodal_indx_end_onring,1);
        Node_e_z(nodal_indx_start_onring:nodal_indx_end_onring,1) = Node(nodal_indx_start_onring:nodal_indx_end_onring,5);
        Node_e_y(nodal_indx_start_onring:nodal_indx_end_onring,1) = real(((1-(Node_e_z(nodal_indx_start_onring:nodal_indx_end_onring,1)-z_prime).^2/c_ellipsoid^2-Node_e_x(nodal_indx_start_onring:nodal_indx_end_onring,1).^2/a_ellipsoid^2)*b_ellipsoid^2).^0.5);
    end
end

if n_r > k_cal_cr
    i_nc_now = k_cal_cr+1;
    nodal_indx_start_onring = (1+i_nc_now-1)/2*(i_nc_now-1)*n_s+2;

    N_chord = i_cal_cr*n_s;
    [nodes_e,nodes_s,theta_e] = ell_cal();
    indp_va_rim = ([1:1:N_chord]'-1)/N_chord;
    fun_val_rim = theta_e;
    [fitobject_rim,gof_rim] = fit([indp_va_rim;1],[fun_val_rim;1],'cubicspline');
    theta_rotation = atan2(Node(nodal_indx_start_onring:n_t,4),Node(nodal_indx_start_onring:n_t,3));
    for kkjj = 1:1:length(theta_rotation)
        if theta_rotation(kkjj)+2*pi() < 2*pi()
            theta_rotation(kkjj) = theta_rotation(kkjj)+2*pi();
        end
    end
    theta_ellipse = feval(fitobject_rim,theta_rotation/2/pi())*2*pi();
    Node_e_x(nodal_indx_start_onring:n_t,1) = cos(theta_ellipse).*a_temp(nodal_indx_start_onring:n_t,1);
    Node_e_z(nodal_indx_start_onring:n_t,1) = Node(nodal_indx_start_onring:n_t,5);
    Node_e_y(nodal_indx_start_onring:n_t,1) = real(((1-(Node_e_z(nodal_indx_start_onring:n_t,1)-z_prime).^2/c_ellipsoid^2-Node_e_x(nodal_indx_start_onring:n_t,1).^2/a_ellipsoid^2)*b_ellipsoid^2).^0.5);
end

if n_r < k_cal_cr+1 || n_r < 0.5*k_cal_cr*(c_b+1)+2
    i_nc_now = n_r;
    nodal_indx_start_onring = (1+i_nc_now-1)/2*(i_nc_now-1)*n_s+2;
    
    N_chord = (i_nc_now)*n_s;
    [nodes_e,nodes_s,theta_e] = ell_cal();
    indp_va_rim = ([1:1:N_chord]'-1)/N_chord;
    fun_val_rim = theta_e;
    [fitobject_rim,gof_rim] = fit([indp_va_rim;1],[fun_val_rim;1],'cubicspline');
    theta_rotation = atan2(Node(nodal_indx_start_onring:n_t,4),Node(nodal_indx_start_onring:n_t,3));
    for kkjj = 1:1:length(theta_rotation)
        if theta_rotation(kkjj)+2*pi() < 2*pi()
            theta_rotation(kkjj) = theta_rotation(kkjj)+2*pi();
        end
    end
    theta_ellipse = feval(fitobject_rim,theta_rotation/2/pi())*2*pi();
    Node_e_x(nodal_indx_start_onring:n_t,1) = cos(theta_ellipse).*a_temp(nodal_indx_start_onring:n_t,1);
    Node_e_z(nodal_indx_start_onring:n_t,1) = Node(nodal_indx_start_onring:n_t,5);
    Node_e_y(nodal_indx_start_onring:n_t,1) = real(((1-(Node_e_z(nodal_indx_start_onring:n_t,1)-z_prime).^2/c_ellipsoid^2-Node_e_x(nodal_indx_start_onring:n_t,1).^2/a_ellipsoid^2)*b_ellipsoid^2).^0.5);
end

% figure
% plot(indp_va_rim,fun_val_rim,'.');hold on
% plot(fitobject_rim,'-r');hold off

for ii_sign = 1:1:length(Node_e_y)
    if Node_e_y(ii_sign)*Node(ii_sign,4) < 0 
        Node_e_y(ii_sign) = - Node_e_y(ii_sign);
    end
end
Node_e = [Node_e_x Node_e_y Node_e_z];
%% Nodal limitation
if length(Node_e(:,1)) > 200
    error('The Toolbox is limited up to 200 nodes in one reflector surface.')
end
%% Step 6 Projections
Node_design = zeros(n_t,3);
if c_b == 0
    n_nobc = n_t-n_r*n_s;
else
    n_nobc = n_t-(n_r-1)*n_s/((c_b+1)/2);
end   
syms x
kk = 1;
xi_p = Node_e_x;
yi_p = Node_e_y;
zi_p = Node_e_z;

zp_ini = -1.8;
zp_step = 0.2;
zp_end = 3;
temp_zp = 1+length([zp_ini:zp_step:zp_end]);
for zp_span = zp_ini:zp_step:zp_end
    z_p(temp_zp-kk) = -10^zp_span*F;

    b2 = 1/4/F_parent*((cos(angle_off)*xi_p+sin(angle_off)*(zi_p-z_p(temp_zp-kk))).^2+yi_p.^2);
    b1 = -sin(angle_off)*xi_p+cos(angle_off)*(zi_p-z_p(temp_zp-kk))+1/2/F_parent*(cos(angle_off)*xi_p+sin(angle_off)*(zi_p-z_p(temp_zp-kk)))*(sin(angle_off)*z_p(temp_zp-kk)+x_o);
    b0 = cos(angle_off)*z_p(temp_zp-kk)+z_o+1/4/F_parent*(sin(angle_off)*z_p(temp_zp-kk)+x_o)^2;


    yita = zeros(n_t,1);
    sol =[(-b1 + (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2), (-b1 - (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2)];
    temp_z = sol(:,1).*(zi_p-z_p(temp_zp-kk))+z_p(temp_zp-kk);
    temp_yita = sol(:,1);
    for i = 1:1:n_t
        if temp_z(i)>-1e-8 && temp_yita(i)>-1e-8
            yita(i) = sol(i,1);
        else
            yita(i) = sol(i,2);
        end
    end
    Node_design = [yita.*xi_p,yita.*yi_p,max((yita.*(zi_p-z_p(temp_zp-kk))+z_p(temp_zp-kk)),0)];
    if min(abs(b2)) <1e-8
        [temptemp, A_indx] = min(abs(b2));
        yita_special = -b0/b1(A_indx);
        Node_design(A_indx,:) = [xi_p(A_indx)*yita_special yi_p(A_indx)*yita_special (zi_p(A_indx)-z_p(kk))*yita_special+z_p(kk)];
    end

    MemberL = zeros(length(M_indx(:,1)),1);
    for iii = 1: length(M_indx(:,1))
        i0 = M_indx(iii,1);
        iL = M_indx(iii, 2);
        dx = Node_design(iL, 1) - Node_design(i0, 1);
        dy = Node_design(iL, 2) - Node_design(i0, 2);
        dz = Node_design(iL, 3) - Node_design(i0, 3);
        MemberL(iii) = sqrt(dx^2 + dy^2 + dz^2);
        if abs(MemberL(iii)) < 1e-6, 
            error('Warning: the length of a member = 0!'); 
        end
    end
    L_tt(temp_zp-kk) = sum(MemberL);
    error_mrms_nobc_p(temp_zp-kk) = 1/16/15^0.5*max(MemberL(1:m_nobc))^2/F;
    error_mrms_p(temp_zp-kk) = 1/16/15^0.5*max(MemberL)^2/F;
    uniform_ratio_p(temp_zp-kk) = max(MemberL)/min(MemberL);
    uniform_ratio_nobc_p(temp_zp-kk) = max(MemberL(1:m_nobc))/min(MemberL(1:m_nobc));
    temp_a = (Node_design(n_nobc-n_s*(n_r-1)+1,1)^2+Node_design(n_nobc-n_s*(n_r-1)+1,2)^2)^0.5;
    Deff_tt(temp_zp-kk) = (temp_a*cos(angle_off)+abs(x_o))*2;
    kk = kk + 1;
end

for zp_span = 0:0.05:1
    z_p(kk) = zp_span*0.9*H_p;
    xi_p = Node_e_x;
    yi_p = Node_e_y;
    zi_p = Node_e_z;
    b2 = 1/4/F_parent*((cos(angle_off)*xi_p+sin(angle_off)*(zi_p-z_p(kk))).^2+yi_p.^2);
    b1 = -sin(angle_off)*xi_p+cos(angle_off)*(zi_p-z_p(kk))+1/2/F_parent*(cos(angle_off)*xi_p+sin(angle_off)*(zi_p-z_p(kk)))*(sin(angle_off)*z_p(kk)+x_o);
    b0 = cos(angle_off)*z_p(kk)+z_o+1/4/F_parent*(sin(angle_off)*z_p(kk)+x_o)^2;

    yita = zeros(n_t,1);
    sol =[(-b1 + (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2), (-b1 - (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2)];
    temp_z = sol(:,1).*(zi_p-z_p(kk))+z_p(kk);
    temp_yita = sol(:,1);
    for i = 1:1:n_t
        if temp_z(i)>-1e-8 && temp_yita(i)>-1e-8
            yita(i) = sol(i,1);
        else
            yita(i) = sol(i,2);
        end
    end
    Node_design = [yita.*xi_p,yita.*yi_p,max((yita.*(zi_p-z_p(kk))+z_p(kk)),0)];
    if min(abs(b2)) <1e-8
        [temptemp, A_indx] = min(abs(b2));
        yita_special = -b0/b1(A_indx);
        Node_design(A_indx,:) = [xi_p(A_indx)*yita_special yi_p(A_indx)*yita_special (zi_p(A_indx)-z_p(kk))*yita_special+z_p(kk)];
    end
    MemberL = zeros(length(M_indx(:,1)),1);
    for iii = 1: length(M_indx(:,1))
        i0 = M_indx(iii,1);
        iL = M_indx(iii, 2);
        dx = Node_design(iL, 1) - Node_design(i0, 1);
        dy = Node_design(iL, 2) - Node_design(i0, 2);
        dz = Node_design(iL, 3) - Node_design(i0, 3);
        MemberL(iii) = sqrt(dx^2 + dy^2 + dz^2);
        if abs(MemberL(iii)) < 1e-6, 
            error('Warning: the length of a member = 0!'); 
        end
    end
    L_tt(kk) = sum(MemberL);
    error_mrms_nobc_p(kk) = 1/16/15^0.5*max(MemberL(1:m_nobc))^2/F;
    error_mrms_p(kk) = 1/16/15^0.5*max(MemberL)^2/F;
    uniform_ratio_p(kk) = max(MemberL)/min(MemberL);
    uniform_ratio_nobc_p(kk) = max(MemberL(1:m_nobc))/min(MemberL(1:m_nobc));
    temp_a = (Node_design(n_nobc-n_s*(n_r-1)+1,1)^2+Node_design(n_nobc-n_s*(n_r-1)+1,2)^2)^0.5;
    Deff_tt(kk) = (temp_a*cos(angle_off)+abs(x_o))*2;
    kk = kk + 1;
end

[ttttmmmm, z_p_indx] = min(error_mrms_p);
L_t_min = L_tt(z_p_indx);

if flag_effect == 0
    z_p_final = z_p(z_p_indx);

    xi_p = Node_e_x;
    yi_p = Node_e_y;
    zi_p = Node_e_z;
    b2 = 1/4/F_parent*((cos(angle_off)*xi_p+sin(angle_off)*(zi_p-z_p_final)).^2+yi_p.^2);
    b1 = -sin(angle_off)*xi_p+cos(angle_off)*(zi_p-z_p_final)+1/2/F_parent*(cos(angle_off)*xi_p+sin(angle_off)*(zi_p-z_p_final))*(sin(angle_off)*z_p_final+x_o);
    b0 = cos(angle_off)*z_p_final+z_o+1/4/F_parent*(sin(angle_off)*z_p_final+x_o)^2;

    yita = zeros(n_t,1);
    sol =[(-b1 + (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2), (-b1 - (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2)];
    temp_z = sol(:,1).*(zi_p-z_p_final)+z_p_final;
    temp_yita = sol(:,1);
    for i = 1:1:n_t
        if temp_z(i)>-1e-8 && temp_yita(i)>-1e-8
            yita(i) = sol(i,1);
        else
            yita(i) = sol(i,2);
        end
    end
    Node_design = [yita.*xi_p,yita.*yi_p,max((yita.*(zi_p-z_p_final)+z_p_final),0)];
    if min(abs(b2)) <1e-8
        [temptemp, A_indx] = min(abs(b2));
        yita_special = -b0/b1(A_indx);
        Node_design(A_indx,:) = [xi_p(A_indx)*yita_special yi_p(A_indx)*yita_special (zi_p(A_indx)-z_p_final)*yita_special+z_p_final];
    end
else
%         figure
%         plotyy(error_mrms_nobc_p,L_tt,error_mrms_nobc_p,Deff_tt)
    z_p_indx = 0;
    z_p_final = -inf;
    b2 = 1/4/F_parent*sin(angle_off)^2;
    b1= 1/2/F_parent*(cos(angle_off)*xi_p+x_o)*sin(angle_off)+cos(angle_off);
    b0 = 1/4/F_parent*((cos(angle_off)*xi_p+x_o).^2+yi_p.^2)-sin(angle_off)*xi_p+z_o;
    Node_design = [Node_e_x Node_e_y (-b1+(b1.^2-4*b2.*b0).^0.5)./(2*b2)];
end

Node_design_global = [Node_design(:,1).*cos(angle_off)+Node_design(:,3).*sin(angle_off)+x_o Node_design(:,2) -Node_design(:,1).*sin(angle_off)+Node_design(:,3).*cos(angle_off)+z_o];

%% Final Calculations
MemberL_prop = zeros(length(M_indx(:,1)),1);
for kkk = 1: length(M_indx(:,1))
    i0 = M_indx(kkk,1);
    iL = M_indx(kkk, 2);
    dx = Node_design(iL, 1) - Node_design(i0, 1);
    dy = Node_design(iL, 2) - Node_design(i0, 2);
    dz = Node_design(iL, 3) - Node_design(i0, 3);
    MemberL_prop(kkk) = sqrt(dx^2 + dy^2 + dz^2);
    if MemberL_prop(kkk) == 0, 
        error('Warning: the length of a member = 0!'); 
    end
end
L_t = sum(MemberL_prop);
L_t_nobc = sum(MemberL_prop(1:m_nobc));
