function [cut] = getcut(A, idxlist)
  cut=nnz(A);
  nparts = size(idxlist,2);
  for i=1:nparts
    range = idxlist(i,1):idxlist(i,2);
    cut = cut - full(nnz(A(range,range)));
  end
end
