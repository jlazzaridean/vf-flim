function [aFs,tFs] = sortATs(aF,tF)
%utility function written by Julia Lazzari-Dean to sort the outputs of
%Matlab fitting so that the lifetimes are in ascending order.

%transposing to set up an easily sortable array
sortable = zeros(size(aF,2),2);
sortable(:,1) = tF';
sortable(:,2) = aF';

%sort this combined table
sorted = sortrows(sortable);

%read out the values, still paired and now in order
aFs = sorted(:,2)';
tFs = sorted(:,1)';

end

