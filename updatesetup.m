%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updatesetup(F_update, D_update, e_update)

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



    if flag_shape == 1 || flag_shape == 2
        F_ref = F_update;
        D_ref = D_update;
    end
    if flag_shape == 3
        D_parent = D_update;
        F_parent = F_update;
        e_offset = e_update;
    end
    
    flag_update = 1;
