%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear
clc
 
% Setup the desired shape and required operating frequency
F_required = 8;
D_required = 30;
e_required = 5;
freq_op_required = 3*1e+9;
flag_effect = 1;
flag_shape = 3;
rms_budget = 50;
% Initial problem setup
error_rms = initialsetup(F_required, D_required, e_required, freq_op_required, flag_effect, flag_shape,rms_budget);
% Preselection of design parameters
n_s = 4;
c_b = 3;% must be 0 when flag_effect = 0
tao_r = 1.1;
tao_c = 0.1;
rou = 2.5;
% Design attempt
flag_wb = 1;
[n_r, n_c, ratio_zeta, w_b, MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F_p, D_p, e_p] = designattempt(tao_r, n_s, tao_c, rou, c_b, flag_wb, []);
% RMS error evaluation
[error_rms_bestfit, deta_bestfit, F_bestfit] = actual_rms_evalu(n_r, n_s, w_b, Node_design_global, M_indx, MemberL_prop(1));
% Best-fit compensation
[temp1, temp3, temp3, D_ca_bestfit, e_off_bestfit, D_p_bestfit, F_p_mod, D_p_mod, e_off_mod] = bestfitcompe(deta_bestfit,F_bestfit,F_p,D_p,e_p);
updatesetup(F_p_mod, D_p_mod, e_off_mod);
% Achieve the final design
[MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F_p_mod, D_p_mod, e_p_mod] = AutoMesh_offset(n_r, n_s, n_c, rou, c_b, ratio_zeta, w_b);
% Measure the actual surface RMS error
[error_rms_mod, deta_bestfit_mod, F_bestfit_mod] = actual_rms_evalu(n_r, n_s, w_b, Node_design_global, M_indx,MemberL_prop(1));
% Evaluate the shape of best-fit surface
[temp1, temp3, temp3, D_ca_bestfit_mod, e_off_bestfit_mod, D_p_bestfit_mod] = bestfitcompe(deta_bestfit_mod,F_bestfit_mod,F_p_mod,D_p_mod,e_p_mod);
    
% Setup the reflector structure for plotting
setup_struc_parameter(M_indx,Node_design, B_C);
% Plot the desired profile
figure
plot_desired(0)
axis equal
view([0 0 1])
% The end of the toolbox demo
disp 'The toolbox works properly.'
