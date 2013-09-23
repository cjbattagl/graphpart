function [B] = partbalanced(A,p,ord)
  n = size(A,1);

  if (nargin==2)
    ord = 1:n;
  else
    assert(length(ord)==n);
  end

  perm = [];
  for i=1:p
    neword = ord(i:(p-1):n);
    perm = [perm neword];
  end
  B=A(perm,perm);
end
