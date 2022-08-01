%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2013 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function memindx_facet = finmemindex_facet(mem_indx_first,M_indx)
clear memindx_facet
i_mem_indx_facet = 1;
node_first = M_indx(mem_indx_first,1);
node_second = M_indx(mem_indx_first,2);
kk_ff = 1;
kk_ff_2= 1;
for i_f_check = 1:1:length(M_indx(:,1))
    if min(abs(node_first-M_indx(i_f_check,:)))< 1e-6
        first_temp(kk_ff)= i_f_check;
        kk_ff = kk_ff+1;
    end
    if min(abs(node_second-M_indx(i_f_check,:)))< 1e-6
        second_temp(kk_ff_2)= i_f_check;
        kk_ff_2 = kk_ff_2+1;
    end
end

third_temp_1 =  (abs(M_indx(first_temp,1) - node_first) < 1e-6).*M_indx(first_temp,2)+(1-(abs(M_indx(first_temp,1) - node_first) < 1e-6)).*M_indx(first_temp,1);
third_temp_2 =  (abs(M_indx(second_temp,1) - node_second) < 1e-6).*M_indx(second_temp,2)+(1-(abs(M_indx(second_temp,1) - node_second) < 1e-6)).*M_indx(second_temp,1);

for i_first = 1:1:length(first_temp)
    for i_second = 1:1:length(second_temp)
%         clear third_temp_1
%         clear third_temp_2
%         if abs(M_indx(first_temp(i_first),1) - node_first) < 1e-6
%             third_temp_1 = M_indx(first_temp(i_first),2);
%         else
%             third_temp_1 = M_indx(first_temp(i_first),1);
%         end
%         if abs(M_indx(second_temp(i_second),1) - node_second) < 1e-6
%             third_temp_2 = M_indx(second_temp(i_second),2);
%         else
%             third_temp_2 = M_indx(second_temp(i_second),1);
%         end        
        if abs(third_temp_1(i_first) - third_temp_2(i_second)) < 1e-6
            memindx_facet(i_mem_indx_facet,:) = [mem_indx_first first_temp(i_first) second_temp(i_second)];
            i_mem_indx_facet = i_mem_indx_facet +1;
        end
    end
end
clear first_temp
clear second_temp


