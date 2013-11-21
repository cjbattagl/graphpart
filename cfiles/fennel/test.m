function [B A cut] = test(fname)
  a=load('parts.mat');
  A=load(fname);
  A=A.Problem.A;
  nparts = max(a(:,2));
  P = [];
  idxlist = zeros(nparts,2); %contains indices of partition matrix
  psizes = zeros(nparts,1);
  N = size(A,1);
  assert(N==size(a,1));

  for i=1:nparts
    [t idx] = find(a(:,2)==i);
    lo = length(P) + 1;
    P = [P; t];
    hi = lo + length(t) - 1;
    idxlist(i,:) = [lo hi];
    psizes(i) = length(t);
  end
  assert(length(P)==N);
  B=A(P,P);
  assert(nnz(A)==nnz(B));
  cut=nnz(B);
  for i=1:nparts
    range = idxlist(i,1):idxlist(i,2);
    cut = cut - nnz(B(range,range));
  end

end
