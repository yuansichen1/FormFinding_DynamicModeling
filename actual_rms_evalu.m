%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2023 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 3.0 $  $Date: 2023/07/07 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [error_rms, deta_final, F_bestfit] = actual_rms_evalu(n_r, n_s, w_b, Node_inuse, M_indx_full, center_lengthmem)

global F_ref
global D_ref
global F_parent
global flag_effect
global flag_shape

if length(Node_inuse(:,1)) > 200
    error('The Toolbox is limited up to 200 nodes in one reflector surface.')
end
if flag_shape == 1
    z_o_act = 0;
    F_focal = F_ref;
else if flag_shape == 2
        z_o_act = D_ref^2/16/F_ref;
        F_focal = F_ref;
    else
        z_o_act = 0;
        F_focal = F_parent;
    end
end
        
% M_indx here is the member ends of the members of working surface, not the full surface
    global F_focal_orginal
    global B_co_ve
    global C_co_ve
    global z_o_center
    global z_c_ve
    global r_c_ve
    global f_ve
    global R_eta
    global R_zeta
    global S_final
    
if flag_effect == 0
    n_cosider = n_r;
    M_indx = M_indx_full;
else
    n_cosider = n_r-1;
    M_indx = M_indx_full(1:(3*n_cosider+1)/2*n_cosider*n_s,:);
end
%%
facet_nodal = finnodeindex_facet(M_indx);
if abs(length(facet_nodal(:,1))-n_cosider^2*n_s*3) > 1e-6
    error('facet index is wrong')
end
Node_first = facet_nodal(:,1);
Node_inuseecond = facet_nodal(:,2);
Node_third = facet_nodal(:,3);
%%

a_co_ve = (Node_inuse(Node_inuseecond,2)-Node_inuse(Node_first,2)).*(Node_inuse(Node_third,3)-Node_inuse(Node_first,3))-(Node_inuse(Node_third,2)-Node_inuse(Node_first,2)).*(Node_inuse(Node_inuseecond,3)-Node_inuse(Node_first,3));
b_co_ve = (Node_inuse(Node_inuseecond,3)-Node_inuse(Node_first,3)).*(Node_inuse(Node_third,1)-Node_inuse(Node_first,1))-(Node_inuse(Node_third,3)-Node_inuse(Node_first,3)).*(Node_inuse(Node_inuseecond,1)-Node_inuse(Node_first,1));
c_co_ve = (Node_inuse(Node_inuseecond,1)-Node_inuse(Node_first,1)).*(Node_inuse(Node_third,2)-Node_inuse(Node_first,2))-(Node_inuse(Node_third,1)-Node_inuse(Node_first,1)).*(Node_inuse(Node_inuseecond,2)-Node_inuse(Node_first,2));
% d_co_ve = -(a_co_ve.*Node_inuse(Node_first,1)+b_co_ve.*Node_inuse(Node_first,2)+c_co_ve.*Node_inuse(Node_first,3));
C_co_ve = -a_co_ve./c_co_ve;
B_co_ve = -b_co_ve./c_co_ve;
% A_co_ve = -d_co_ve./c_co_ve;
z_c_ve = 1/3*(Node_inuse(Node_first,3)+Node_inuse(Node_inuseecond,3)+Node_inuse(Node_third,3));
y_c_ve = 1/3*(Node_inuse(Node_first,2)+Node_inuse(Node_inuseecond,2)+Node_inuse(Node_third,2));
x_c_ve = 1/3*(Node_inuse(Node_first,1)+Node_inuse(Node_inuseecond,1)+Node_inuse(Node_third,1));
%%
% sphere
if flag_shape == 1
    R_eta = ones(n_cosider^2*n_s*3,1)*2*F_focal;
    R_zeta = ones(n_cosider^2*n_s*3,1)*2*F_focal;
else
% parabolic
    r_c_ve = (y_c_ve.^2+x_c_ve.^2).^0.5;
    R_eta = 2*F_focal*(1+(r_c_ve/2/F_focal).^2).^(3/2);
    R_zeta = r_c_ve.*(1+(2*F_focal./r_c_ve).^2).^0.5;
end
%%
gama_ve = acos((B_co_ve.*y_c_ve+C_co_ve.*x_c_ve)./((B_co_ve.*y_c_ve+C_co_ve.*x_c_ve).^2+y_c_ve.^2+x_c_ve.^2).^0.5);
theta_ve = atan2(x_c_ve,y_c_ve);
eta_2_local_ve = ((Node_inuse(Node_inuseecond,2)-Node_inuse(Node_first,2)).*cos(theta_ve)+(Node_inuse(Node_inuseecond,1)-Node_inuse(Node_first,1)).*sin(theta_ve))./sin(gama_ve);
eta_3_local_ve = ((Node_inuse(Node_third,2)-Node_inuse(Node_first,2)).*cos(theta_ve)+(Node_inuse(Node_third,1)-Node_inuse(Node_first,1)).*sin(theta_ve))./sin(gama_ve);
l_2_length = ((Node_inuse(Node_inuseecond,1)-Node_inuse(Node_first,1)).^2+(Node_inuse(Node_inuseecond,2)-Node_inuse(Node_first,2)).^2+(Node_inuse(Node_inuseecond,3)-Node_inuse(Node_first,3)).^2).^0.5;
l_3_length = ((Node_inuse(Node_third,1)-Node_inuse(Node_first,1)).^2+(Node_inuse(Node_third,2)-Node_inuse(Node_first,2)).^2+(Node_inuse(Node_third,3)-Node_inuse(Node_first,3)).^2).^0.5;
l_23_length = ((Node_inuse(Node_third,1)-Node_inuse(Node_inuseecond,1)).^2+(Node_inuse(Node_third,2)-Node_inuse(Node_inuseecond,2)).^2+(Node_inuse(Node_third,3)-Node_inuse(Node_inuseecond,3)).^2).^0.5;
yeta_2 = (abs(l_2_length.^2-eta_2_local_ve.^2)).^0.5;
yeta_3_temp = -(abs((l_3_length.^2-eta_3_local_ve.^2))).^0.5;
S_temp = abs(eta_2_local_ve.*yeta_3_temp-eta_3_local_ve.*yeta_2)/2;
length_aveg = (l_2_length+l_3_length+l_23_length)/2;
S_temp_2 = (length_aveg.*(length_aveg-l_2_length).*(length_aveg-l_3_length).*(length_aveg-l_23_length)).^0.5;
yeta_3 = (-1).^(abs(S_temp_2-S_temp)./S_temp > 1e-3).*yeta_3_temp;
S_final = abs(eta_2_local_ve.*yeta_3-eta_3_local_ve.*yeta_2)/2;
f_ve = (eta_2_local_ve.^2-eta_2_local_ve.*eta_3_local_ve+eta_3_local_ve.^2)./R_eta+(yeta_2.^2-yeta_2.*yeta_3+yeta_3.^2)./R_zeta;
if max(abs(imag(f_ve))) > 1e-12
    error('f_ve has imaginary part, program stopped')
else
    f_ve = real(f_ve);
end
%%
if flag_shape == 1
    % sphere
    if max(abs(imag(f_ve)))==0 && max(abs(imag(S_final)))==0 && max(abs(imag(R_eta)))==0 && max(abs(imag(R_zeta)))==0
        a_qfunc = real(sum(S_final));
        b_qfunc = real(sum(-S_final.*f_ve/6));
        c_qfunc = real(sum(S_final.*(f_ve.^2/120-S_final.^2./(90.*R_eta.*R_zeta))));
        error_rms = real(((4*a_qfunc*c_qfunc-b_qfunc^2)/4/a_qfunc/sum(S_final))^0.5);
        deta_final = real(-b_qfunc./2./a_qfunc);
        F_bestfit = F_focal;
    else
        error('imag problem, program stopped')
    end
    % error_phi_ve = S_final.*(deta_temp.^2-deta_temp.*f_ve/6+f_ve.^2/120-S_final.^2./(90.*R_eta.*R_zeta));
    % error_rms = (sum(error_phi_ve)/sum(S_final))^0.5;
else
    % parabolic
    H_maxmax = 2*F_focal-((2*F_focal)^2-center_lengthmem^2/3)^0.5;
    
    z_o_center = z_o_act;
    F_focal_orginal = F_focal;
    [outppp,error_rms] = fminunc(@rms_parabolic_search,[F_focal,3/4*H_maxmax]);
    F_bestfit = outppp(1);
    deta_final = outppp(2);
end