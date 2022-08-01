clear
clc
 
% Setup the desired shape and required operating frequency
F_required = 8;
D_required = 10;
e_required = [];
freq_op_required = 2*1e+9;
flag_effect = 1;
flag_shape = 2;
rms_budget = 50;
% Initial problem setup
error_rms = initialsetup(F_required, D_required, e_required, freq_op_required, flag_effect, flag_shape,rms_budget);
% Preselection of design parameters
n_s = 6;
c_b = 5;% must be 0 when flag_effect = 0
tao_r = 1.45;
tao_c = 0.65;
rou = 2.4;
% Design attempt
flag_wb = 1;
[n_r, n_c, ratio_zeta, w_b, MemberL_prop,M_indx, Node_design,Node_design_global, B_C, Indx_node_load, L_t, L_t_nobc, F, D] = designattempt(tao_r, n_s, tao_c, rou, c_b, flag_wb, []);
% RMS error evaluation
[error_rms_bestfit, deta_bestfit, F_bestfit] = actual_rms_evalu(n_r, n_s, w_b, Node_design, M_indx, MemberL_prop(1));
% Best-fit compensation
[D_bestfit, D_mod, F_mod] = bestfitcompe(deta_bestfit,F_bestfit,F,D);
updatesetup(F_mod, D_mod, []);
% Achieve the final design
[MemberL_prop,M_indx, Node_design, B_C, Indx_node_load, L_t, L_t_nobc, F_mod, D_mod] = AutoMesh(n_r, n_s, n_c, rou, c_b, ratio_zeta, w_b);
% Measure the actual surface RMS error
[error_rms_mod, deta_bestfit_mod, F_bestfit_mod] = actual_rms_evalu(n_r, n_s, w_b, Node_design, M_indx,MemberL_prop(1));
% Evaluate the shape of best-fit surface
[D_bestfit_mod] = bestfitcompe(deta_bestfit_mod,F_bestfit_mod,F_mod,D_mod,[]);
    
% Setup the reflector structure for plotting
setup_struc_parameter(M_indx,Node_design, B_C);
% Plot the desired profile
figure
plot_desired(0)
 axis equal
view([0 0 1])
% The end of the toolbox demo
disp 'The toolbox works properly.'
