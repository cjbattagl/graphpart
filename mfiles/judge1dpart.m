% Report metrics on 1d partition quality into p equal parts
% future: allow bins to represent partitions
function [balance cut loads] = judge1dpart(A,p)
  n = size(A,1);
  m = nnz(A);
  loads = zeros(p,1);

  for i=1:p
    loads(i) = sum(sum(A((i-1)*floor(n/p)+1:i*floor(n/p),:)));
  end
  balance = max(loads)/min(loads)
  cut = cutsize(A,p)/m;
end
