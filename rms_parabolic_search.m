%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function error_rms = rms_parabolic_search(input)

F_focal_temp=input(1);
deta_temp=input(2);
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

r_q_ve = -2*F_focal_orginal./(B_co_ve.^2+C_co_ve.^2).^0.5+2.*(F_focal_orginal^2./(B_co_ve.^2+C_co_ve.^2)+F_focal_orginal.*(z_o_center-z_c_ve+r_c_ve./((B_co_ve.^2+C_co_ve.^2).^0.5))).^0.5;
r_p_ve = -2*F_focal_temp./(B_co_ve.^2+C_co_ve.^2).^0.5+2.*(F_focal_temp^2./(B_co_ve.^2+C_co_ve.^2)+F_focal_temp.*(z_o_center-z_c_ve-deta_temp+r_c_ve./((B_co_ve.^2+C_co_ve.^2).^0.5))).^0.5;
z_q_ve = z_o_center-r_q_ve.^2/4/F_focal_orginal;
z_p_ve = z_o_center-r_p_ve.^2/4/F_focal_temp;
a_span = ((z_p_ve-z_q_ve).^2+(r_p_ve-r_q_ve).^2).^0.5;
if max(abs(imag(f_ve)))==0 && max(abs(imag(S_final)))==0 && max(abs(imag(R_eta)))==0 && max(abs(imag(R_zeta)))==0
    error_phi_ve = real(S_final.*(a_span.^2-a_span.*f_ve/6+f_ve.^2/120-S_final.^2./(90.*R_eta.*R_zeta)));
else
    error('imag problem, program stopped')
end
error_rms = (sum(error_phi_ve)/sum(S_final))^0.5;