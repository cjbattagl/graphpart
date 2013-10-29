% This function displays the 'cut-size' of the matrix
% if it were to be logically partitioned into k equally-sized partitions
% in the order that the matrix is presented in
% it also displays the percentage of edges in the edge-cut

function f = cutsize(A, k)
  m = size(A,1);
  S = sum(sum(A));
  for i=0:k-1
    idxs = (i*floor(m/k)+1:(i+1)*floor(m/k));
    S = S - sum(sum(A(idxs,idxs)));
  end
  disp(S)              % #edges
  disp(100*S/nnz(A))   % %edges cut
  f = 100*S/nnz(A);
end