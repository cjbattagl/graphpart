%How to use:
% ./fennel -v -e 'ca-AstroPh.mtx' 'MM'
% [B A cut]=test('ca-AstroPh.mat');

function [B A C cut] = test(fname)
  a=load('parts.mat');
  A=load(fname);
  A=spconvert(A);
  maxlen = max(size(A))+1
  A(maxlen,maxlen)=1;
  %A=A.Problem.A;
  nparts = max(a(:,2))-1;
  P = [];
  idxlist = zeros(nparts,2); %contains indices of partition matrix
  psizes = zeros(nparts,1);
  N = size(A,1);
  %assert(N==size(a,1));
  rem = 1:N;
  for i=1:nparts
    [t idx] = find(a(:,2)==i);
    verts = a(t,1);
    lo = length(P) + 1;
    rem(verts) = -1;
    P = [P; verts];
    hi = lo + length(t) - 1;
    idxlist(i,:) = [lo hi];
    psizes(i) = length(t);
  end
  t = find(rem==-1);
  rem(t) = [];
  rem = rem';
  C = A(P,P);
  B=A([P; rem],[P; rem]);

  % Sanity checks
  assert(length(unique(P))==length(P));
  %assert(length(P)==N);
  %assert(nnz(A)==nnz(B));
  
  % Verify cut size
  cut=nnz(C);
  for i=1:nparts
    range = idxlist(i,1):idxlist(i,2);
    cut = cut - full(nnz(C(range,range)));
  end

  cut = cut / nnz(C);

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
