%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2023 by Sichen Yuan
%   Created: 2022/08/01
%   $Revision: 2.0 $  $Date: 2023/07/07 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n_r, n_c, zeta_final, wb_final, MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = designattempt(tao_r, n_s, tao_c, rou, c_b, flag_wb, w_b)

% flag_wb = 0 wb is assigned; 1: calculated

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

n_r =round(tao_r*n_r_bar);
n_c = round(tao_c*n_r);
%%
Node_design_global = [];
e_p = [];
if flag_wb == 0
    if flag_effect == 0
        if flag_shape == 1 || flag_shape == 2
            ratio_zeta = 1;
            w_b = 0;
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        else
            ratio_zeta = 1;
            w_b = 0;
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        end
    else
        if flag_shape == 1 || flag_shape == 2
            ratio_zeta = 1;
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        else
            ratio_zeta = 1;
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        end
    end
else
    if flag_effect == 0
        if flag_shape == 1 || flag_shape == 2
            ratio_zeta = 1;
            w_b = 0;
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        else
            ratio_zeta = 1;
            w_b = 0;
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        end
    else
        if flag_shape == 1 || flag_shape == 2
            ratio_zeta = 1;
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
    
            w_b = max(ceil((1/asin(D_ref_eff/4/F_ref_eff)*asin(D_ref_eff/4/F_ref_eff/cos(pi()*(c_b-1)/2/(n_r-1)/n_s))-1)*n_i_i_temp*10)/10,1);
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            w_b = max(ceil((1/asin(D_ref_eff/4/F_ref_eff)*asin(D_ref_eff/4/F_ref_eff/cos(pi()*(c_b-1)/2/(n_r-1)/n_s))-1)*n_i_i_temp*10)/10,1);
            [MemberL_prop, M_indx, Node_design,B_C, Indx_node_load, L_t, L_t_nobc, F, D] = AutoMesh(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        else
            ratio_zeta = 1;
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
    
            w_b = max(ceil((1/asin(D_ref_eff/4/F_ref_eff)*asin(D_ref_eff/4/F_ref_eff/cos(pi()*(c_b-1)/2/(n_r-1)/n_s))-1)*n_i_i_temp*10)/10,1);
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
            ratio_zeta = max(MemberL_prop((3*(n_r-2)+1)/2*(n_r-2)*n_s+1:(3*(n_r-1)+1)/2*(n_r-1)*n_s))./max(MemberL_prop(1:12));
            w_b = max(ceil((1/asin(D_ref_eff/4/F_ref_eff)*asin(D_ref_eff/4/F_ref_eff/cos(pi()*(c_b-1)/2/(n_r-1)/n_s))-1)*n_i_i_temp*10)/10,1);
            [MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D, e_p] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b,ratio_zeta,w_b);
        end
    end
end

zeta_final = ratio_zeta;
wb_final = w_b;
%% Nodal limitation
if length(Node_design(:,1)) > 200
    error('The Toolbox is limited up to 200 nodes in one reflector surface.')
end
       