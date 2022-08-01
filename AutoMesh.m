%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b, ratio_zeta, w_b)

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
        F = F_ref;
        D = D_ref;
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
        D = D_ref;
        F_ref = F_ref_eff;
        F = F_ref;
    end
else
    F = F_ref;
    D = D_ref;    
end
%% Step 1 Calculate the geometry and Set up the reference sphere
if flag_shape == 1
    R_s = 2*F;
    H = R_s-(R_s^2-D^2/4)^0.5;
else
    a = 1/4/F;
    H_p = D^2/16/F;
    R_s = 2*F;
    H = R_s-(R_s^2-D^2/4)^0.5;
end
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
%% Nodal limitation
if length(Node(:,1)) > 200
    error('The Toolbox is limited up to 200 nodes in one reflector surface.')
end
%% Step 6 Projections
Node_design = zeros(n_t,3);
if flag_shape == 2
    syms x
    kk = 1;
    xi_p = Node(:,3);
    yi_p = Node(:,4);
    zi_p = Node(:,5);
%     Node_design = [xi_p, yi_p, H_p-a*(xi_p.^2+yi_p.^2)];
%     for iii = 1: length(M_indx(:,1))
%         i0 = M_indx(iii,1);
%         iL = M_indx(iii, 2);
%         dx = Node_design(iL, 1) - Node_design(i0, 1);
%         dy = Node_design(iL, 2) - Node_design(i0, 2);
%         dz = Node_design(iL, 3) - Node_design(i0, 3);
%         MemberL(iii) = sqrt(dx^2 + dy^2 + dz^2);
%         if abs(MemberL(iii)) < 1e-6, 
%             error('Warning: the length of a member = 0!'); 
%         end
%     end    
%     z_p(kk) = -inf;
%     L_tt(kk) = sum(MemberL);
%     error_mrms_nobc_p(kk) = 1/16/15^0.5*max(MemberL(1:m_nobc))^2/F;
%     error_mrms_p(kk) = 1/16/15^0.5*max(MemberL)^2/F;    
    
    zp_ini = -1.8;
    zp_step = 0.2;
    zp_end = 3;
    temp_zp = 1+length([zp_ini:zp_step:zp_end]);
    for zp_span = zp_ini:zp_step:zp_end
        z_p(temp_zp-kk) = -10^zp_span*F;
        
        xi_p = Node(:,3);
        yi_p = Node(:,4);
        zi_p = Node(:,5);
        b2 = a*(xi_p.^2+yi_p.^2);
        b1 = zi_p-z_p(temp_zp-kk);
        b0 = (z_p(temp_zp-kk)-H_p)*ones(length(Node(:,3)),1);
        yita = zeros(n_t,1);
        sol =[(-b1 + (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2), (-b1 - (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2)];
        temp_z = sol(:,1).*(zi_p-z_p(temp_zp-kk))+z_p(temp_zp-kk);
        temp_yita = sol(:,1);
        for i = 1:1:n_t
            if temp_z(i)>-1e-8 && temp_yita(i)>1e-8
                yita(i) = sol(i,1);
            else
                yita(i) = sol(i,2);
            end
        end
        Node_design = [yita.*xi_p,yita.*yi_p,(yita.*(zi_p-z_p(temp_zp-kk))+z_p(temp_zp-kk))];
        if min(abs(b2)) <1e-8
            [temptemp, A_indx] = min(abs(b2));
            Node_design(A_indx,:) = [0 0 H_p];
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
        kk = kk + 1;
    end
    
    for zp_span = 0:0.05:1
        z_p(kk) = zp_span*0.9*H_p;
        xi_p = Node(:,3);
        yi_p = Node(:,4);
        zi_p = Node(:,5);
        b2 = a*(xi_p.^2+yi_p.^2);
        b1 = zi_p-z_p(kk);
        b0 = (z_p(kk)-H_p)*ones(length(Node(:,3)),1);
        yita = zeros(n_t,1);
        sol =[(-b1 + (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2), (-b1 - (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2)];
        temp_z = sol(:,1).*(zi_p-z_p(kk))+z_p(kk);
        temp_yita = sol(:,1);
        for i = 1:1:n_t
            if temp_z(i)>-1e-8 && temp_yita(i)>1e-8
                yita(i) = sol(i,1);
            else
                yita(i) = sol(i,2);
            end
        end
        Node_design = [yita.*xi_p,yita.*yi_p,(yita.*(zi_p-z_p(kk))+z_p(kk))];
        if min(abs(b2)) <1e-8
            [temptemp, A_indx] = min(abs(b2));
            Node_design(A_indx,:) = [0 0 H_p];
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
        kk = kk + 1;
    end
    
    [ttttmmmm, z_p_indx] = min(error_mrms_p);
    L_t_min = L_tt(z_p_indx);
        
    if flag_effect == 0
        z_p_final = z_p(z_p_indx);
  
        xi_p = Node(:,3);
        yi_p = Node(:,4);
        zi_p = Node(:,5);
        b2 = a*(xi_p.^2+yi_p.^2);
        b1 = zi_p-z_p_final;
        b0 = (z_p_final-H_p)*ones(length(Node(:,3)),1);
        yita = zeros(n_t,1);
        sol =[(-b1 + (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2), (-b1 - (b1.^2 - 4.*b0.*b2).^(1/2))./(2.*b2)];
        temp_z = sol(:,1).*(zi_p-z_p_final)+z_p_final;
        temp_yita = sol(:,1);
        for i = 1:1:n_t
            if temp_z(i)>-1e-8 && temp_yita(i)>1e-8
                yita(i) = sol(i,1);
            else
                yita(i) = sol(i,2);
            end
        end
        Node_design = [yita.*xi_p,yita.*yi_p,(yita.*(zi_p-z_p_final)+z_p_final)];
        if min(abs(b2)) <1e-8
            [temptemp, A_indx] = min(abs(b2));
            Node_design(A_indx,:) = [0 0 H_p];
        end
    else
%         plotyy(error_mrms_nobc_p,L_tt,error_mrms_nobc_p,Deff_tt)
        z_p_indx = 0;
        z_p_final = -inf;
        Node_design = [Node(:,3), Node(:,4), H_p-a*(Node(:,3).^2+Node(:,4).^2)];
    end
else
    Node_design = Node(:,3:5);
    z_p_final = -inf;
    z_p = -inf;
    L_tt = -inf;
    error_mrms_nobc_p = [];
    error_mrms_p = [];
end
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