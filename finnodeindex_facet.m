%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (c) 2022 by Sichen Yuan
%   Created: 2022/05/30
%   $Revision: 1.0 $  $Date: 2022/05/30 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function facet_nodal = finnodeindex_facet(M_indx)
i_node_facet=1;
for i_mm = 1:1:length(M_indx(:,1))
    node_first = M_indx(i_mm,1);
    node_second = M_indx(i_mm,2);

    third_temp_1 =  (M_indx(max((M_indx==node_first),[],2),1) == node_first).*M_indx(max((M_indx==node_first),[],2),2)+(M_indx(max((M_indx==node_first),[],2),2) == node_first).*M_indx(max((M_indx==node_first),[],2),1);
    third_temp_2 =  (M_indx(max((M_indx==node_second),[],2),1) == node_second).*M_indx(max((M_indx==node_second),[],2),2)+(M_indx(max((M_indx==node_second),[],2),2) == node_second).*M_indx(max((M_indx==node_second),[],2),1);

    third_check = intersect(third_temp_1,third_temp_2);
    num_third_check = length(third_check);
    if num_third_check > 0
        facet_nodal(i_node_facet:i_node_facet+num_third_check-1,:) = [ones(num_third_check,1)*node_first ones(num_third_check,1)*node_second third_check];
        i_node_facet = i_node_facet + num_third_check;
    end
end

