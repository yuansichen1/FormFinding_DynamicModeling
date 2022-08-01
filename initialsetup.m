%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function error_rms = initialsetup(F_required, D_required, e_required, freq_op_required, flag_effect1, flag_shape1,rms_budget)

% F_required, D_required, e_required in meters and in parent, freq_op_required in Hz.
% flag_shape: 1: sphere 2: center parabolic 3: offset parabolic
% flag_effect: 0: no effective concern 1: effective region

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

flag_effect = flag_effect1;
flag_shape = flag_shape1;

if flag_effect == 0
    if flag_shape == 1
        F_ref = F_required;
        D_ref = D_required;
        Speed_light = 299792458;
        lamda_wl = Speed_light/freq_op_required;
        error_rms = lamda_wl/rms_budget;
        L_f_bar = 4*15^(1/4)*(error_rms/D_ref*F_ref/D_ref)^0.5*D_ref;
        theta_f = 2*asin(L_f_bar/4/F_ref);
        theta = asin(D_ref/4/F_ref);
        n_r_bar = round(theta/theta_f);
    end
    if flag_shape == 2
        F_ref = F_required;
        D_ref = D_required;
        Speed_light = 299792458;
        lamda_wl = Speed_light/freq_op_required;
        error_rms = lamda_wl/rms_budget;
        L_f_bar = 4*15^(1/4)*(error_rms/D_ref*F_ref/D_ref)^0.5*D_ref;
        theta_f = 2*asin(L_f_bar/4/F_ref);
        theta = asin(D_ref/4/F_ref);
        n_r_bar = round(theta/theta_f);
    end
    if flag_shape == 3
        D_parent = D_required;
        F_parent = F_required;
        FoD_parent = F_parent/D_parent;
        e_offset = e_required;
        ratio_off = e_offset/D_parent*2;
        R_c = 1/2*(D_parent/2-e_offset);
        H_parent = D_parent^2/16/F_parent;
        e_zz = ratio_off^2*D_parent^2/16/F_parent;
        z_o = -1/2*e_zz-1/2*H_parent;
        x_o = e_offset+R_c;
        y_o = 0;
        angle_off = atan((1+ratio_off)*D_parent/8/F_parent);
        a_ellipse = R_c/cos(angle_off);
        b_ellipse = R_c;
        D_circle = 2*R_c;
        F_ref = F_parent;
        D_ref = 2*a_ellipse;
        FoD_ref = F_ref/D_ref;
        FoD_aper = F_ref/D_circle;
        Speed_light = 299792458;
        lamda_wl = Speed_light/freq_op_required;
        error_rms = lamda_wl/rms_budget;
        L_f_bar = 4*15^(1/4)*(error_rms/D_ref*F_ref/D_ref)^0.5*D_ref;
        theta_f = 2*asin(L_f_bar/4/F_ref);
        theta = asin(D_ref/4/F_ref);
        n_r_bar = round(theta/theta_f);
    end
else
    if flag_shape == 1
        F_ref_eff = F_required;
        D_ref_eff = D_required;
        Speed_light = 299792458;
        lamda_wl = Speed_light/freq_op_required;
        error_rms = lamda_wl/rms_budget;
        L_f_bar = 4*15^(1/4)*(error_rms/D_ref_eff*F_ref_eff/D_ref_eff)^0.5*D_ref_eff;
        theta_f = 2*asin(L_f_bar/4/F_ref_eff);
        theta = asin(D_ref_eff/4/F_ref_eff);
        n_r_bar = round(theta/theta_f)+1;
    end
    if flag_shape == 2
        F_ref_eff = F_required;
        D_ref_eff = D_required;
        Speed_light = 299792458;
        lamda_wl = Speed_light/freq_op_required;
        error_rms = lamda_wl/rms_budget;
        L_f_bar = 4*15^(1/4)*(error_rms/D_ref_eff*F_ref_eff/D_ref_eff)^0.5*D_ref_eff;
        theta_f = 2*asin(L_f_bar/4/F_ref_eff);
        theta = asin(D_ref_eff/4/F_ref_eff);
        n_r_bar = round(theta/theta_f)+1;
    end
    if flag_shape == 3
        D_p_eff = D_required;
        F_p_eff = F_required;
        e_eff = e_required;
        ratio_off_eff = e_eff/D_p_eff*2;
        R_c_eff = 1/2*(D_p_eff/2-e_eff);
        H_p_eff = D_p_eff^2/16/F_p_eff;
        e_zz_eff = ratio_off_eff^2*D_p_eff^2/16/F_p_eff;
        z_o_eff = -1/2*e_zz_eff-1/2*H_p_eff;
        x_o_eff = e_eff+R_c_eff;
        y_o_eff = 0;
        angle_eff = atan((1+ratio_off_eff)*D_p_eff/8/F_p_eff);
        a_ellipse_eff = R_c_eff/cos(angle_eff);
        D_circle_eff = 2*R_c_eff;
        D_ref_eff = 2*a_ellipse_eff;
        F_ref_eff = F_p_eff;
        FoD_ref_eff = F_ref_eff/D_ref_eff;
        FoD_aper_eff = F_ref_eff/D_circle_eff;
        Speed_light = 299792458;
        lamda_wl = Speed_light/freq_op_required;
        error_rms = lamda_wl/rms_budget;
        L_f_bar = 4*15^(1/4)*(error_rms/D_ref_eff*F_ref_eff/D_ref_eff)^0.5*D_ref_eff;
        theta_f = 2*asin(L_f_bar/4/F_ref_eff);
        theta = asin(D_ref_eff/4/F_ref_eff);
        n_r_bar = round(theta/theta_f)+1;
    end
end
flag_update = 0;