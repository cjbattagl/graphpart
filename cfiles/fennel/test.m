%How to use:
% ./fennel -v -e 'ca-AstroPh.mtx' 'MM'
% [B A cut]=test('ca-AstroPh.mat');

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

  B=A(P,P);

  % Sanity checks
  assert(length(unique(P))==length(P));
  assert(length(P)==N);
  assert(nnz(A)==nnz(B));
  
  % Verify cut size
  cut=nnz(B);
  for i=1:nparts
    range = idxlist(i,1):idxlist(i,2);
    cut = cut - full(nnz(B(range,range)));
  end

  %%%%%% Other statistics:

  % Compute densities of partitions
  disp('Density of subgraphs:');
  for i=1:nparts
    range = idxlist(i,1):idxlist(i,2);
    dens = full(nnz(B(range,range))) / length(range)^2;
    fprintf('[%d]:\t%1.2f %%\n',i,100*dens);
  end
  disp('');
  % Compute avg degree of partitions
  disp('Average degrees:');
  for i=1:nparts
    range = idxlist(i,1):idxlist(i,2);
    deg = (nnz(B(range,:))/length(range));
    fprintf('[%d]:\t%d\n',i,deg);
  end
  disp('');
end
