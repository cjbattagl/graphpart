% This is a benchmarking experiment where we grow ER(N,p) from t=0 -> infty
% We generate a directed graph with no self-edges
function [lambdas] = lambdaK(N,max_nnz,seed,nparts) 
  iterations = 100; %number of iterations between p_min and p_max
  lbound = N; %default lowerbound on nnz
  bound = 50*10^6; %default upperbound on nnz
  times = []; %vector of benchmark results
  edges = []; %vector of nnz's
  zs = []; %average degrees
  lambdas = [];
  if (nargin > 1)
    if(nargin > 2)
      rand('seed',seed);
    end
    bound = max_nnz;
  end

  p_min = lbound / (N^2 - N); %prob that will gen 'lbound' edges
  p_max = bound / (N^2 - N); %prob that will gen 'bound' edges
  grain = (p_max-p_min)/iterations; 
  x = ones(N,1);
  
  prob = p_min;
  kron_size = ceil(log2(N));
  %ER = erdos_generator(N,prob); 
  kron_list = kronecker_generator(ceil(log2(N)), 1, 0.57, 0.19, 0.19);
  kron_list = kron_list + 1;
  ER = sparse(kron_list(1, :), kron_list(2,:), ones(1,length(kron_list(1,:))), 2^kron_size, 2^kron_size);
  map = metismex('PartGraphKway', ER, nparts); 
  map = (map + 1)';
  map(:,2) = (1:N)';
 
  %P: metis permutation  idxlist: list of indices for each partition 
  [P idxlist] = genperm(map);
  C = ER(P,P);
  cut = getcut(C, idxlist);
  
  assert(nnz(C)==nnz(ER))
  assert(size(P,1)==size(ER,1))

  new_row = [nnz(C) cut cut/nnz(C) nnz(C)/prod(size(C))];
  lambdas = [lambdas; new_row];
  iter = 1;
  for prob = p_min:grain:p_max
    kron_list = kronecker_generator(kron_size, floor(prob*(2^kron_size)), 0.57, 0.19, 0.19);
    kron_list = kron_list + 1;
    ER = sparse(kron_list(1, :), kron_list(2,:), ones(1,length(kron_list(1,:))), 2^kron_size, 2^kron_size);
    if (mod(iter,5) == 0)
       map = metismex('PartGraphKway', ER, nparts);
       map = (map + 1)';
       map(:,2) = (1:N)';
       %P: metis permutation  idxlist: list of indices for each partition 
       [P idxlist] = genperm(map);
    end
    iter = iter + 1;
    ER = ER(P,P);
    cut = getcut(ER,idxlist); 
    new_row = [nnz(ER) cut cut/nnz(ER) nnz(ER)/prod(size(ER))];
    lambdas = [lambdas; new_row];
  end
  disp(num2str(lambdas,'%.2f'))
end
