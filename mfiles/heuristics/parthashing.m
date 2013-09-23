function [B] = parthashing(A,p,ord)
  n = size(A,1);
  if (nargin == 2)
    ord = 1:n;
  else
    assert(length(ord)==n);
  end

  perm = mod(ord, p);
  B=A(perm,perm);
end
