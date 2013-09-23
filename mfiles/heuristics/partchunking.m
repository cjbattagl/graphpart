function [B] = partchunking(A,p,ord)
  if (nargin==2)
    B=A;
  else
    assert(length(ord)=size(A,1));
    B=A(ord,ord);
  end
end
