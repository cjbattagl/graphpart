function [B A cut] = test(fname)
  a=load('parts.mat');
  A=load(fname);
  A=A.Problem.A;
  nparts = max(a(:,2));
  P = [];
  N = size(A,1);
  assert(N==size(a,1));

  for i=1:nparts
    [t idx] = find(a(:,2)==i);
    P = [P; t];
  end
  assert(length(P)==N);
  B=A(P,P);
  cut = 1;
%cutsize

end
