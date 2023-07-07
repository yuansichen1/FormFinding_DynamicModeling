%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2023 by Sichen Yuan
%   Created: 2023/01/01
%   $Revision: 2.0 $  $Date: 2022/07/07 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D_bestfit, D_mod, F_mod, D_ca_bestfit, e_off_bestfit, D_p_bestfit, F_p_mod, D_p_mod, e_off_mod] = bestfitcompe(deta_bestfit,F_bestfit,F_surface,D_surface,e_surface)

%F_surface,D_surface are the parent ones for offset parabolic case

global flag_shape

if flag_shape == 1
    H_surface = 2*F_surface-((2*F_surface)^2-D_surface^2/4)^0.5;
    D_bestfit = ((2*F_bestfit)^2-(2*F_bestfit-H_surface+deta_bestfit)^2)^0.5*2;
    D_mod = D_surface/D_bestfit*D_surface;
    F_mod = F_surface/F_bestfit*F_surface;
    D_ca_bestfit = []; e_off_bestfit =[]; D_p_bestfit=[]; F_p_mod=[]; D_p_mod=[]; e_off_mod=[];
end
if flag_shape == 2
    Hpp = D_surface^2/16/F_surface;
    D_bestfit = (F_bestfit*16*(Hpp-deta_bestfit))^0.5;
    D_mod = D_surface/D_bestfit*D_surface;
    F_mod = F_surface/F_bestfit*F_surface;
    D_ca_bestfit = []; e_off_bestfit =[]; D_p_bestfit=[]; F_p_mod=[]; D_p_mod=[]; e_off_mod=[];
end
if flag_shape == 3
    H_parent = D_surface^2/16/F_surface;
    D_p_prime = (F_bestfit*16*(H_parent-deta_bestfit))^0.5;
    tan_phi_prime = D_p_prime/4/F_bestfit;
    ratio_off = e_surface/D_surface*2;
    tan_phi_off = (1+ratio_off)*D_surface/8/F_surface;
    D_ca = 1/2*D_surface*(1-ratio_off);
    D_ca_bestfit = D_ca-tan_phi_prime/(tan_phi_prime-tan_phi_off)*(D_surface-D_p_prime);
    e_off_bestfit = 1/2*D_surface*ratio_off+1/2*tan_phi_prime/(tan_phi_prime-tan_phi_off)*(D_surface-D_p_prime);
    D_p_bestfit = 1/(tan_phi_prime-tan_phi_off)*(tan_phi_prime*D_p_prime-tan_phi_off*D_surface);
    D_ca_mod = D_ca^2/D_ca_bestfit;
    F_p_mod = F_surface/F_bestfit*F_surface;
    D_p_mod = D_ca_mod-D_ca+D_surface;
    e_off_mod  = 1/2*D_surface*ratio_off-1/2*(D_ca_mod-D_ca);
    D_bestfit=[]; D_mod=[]; F_mod=[];
end
    