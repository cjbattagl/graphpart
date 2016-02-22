% This is a benchmarking experiment where we grow ER(N,p) from t=0 -> infty
% We generate a directed graph with no self-edges
function [lambdas] = lambda1(N,max_nnz,seed,nparts) 
  iterations = 20; %number of iterations between p_min and p_max
  trials = 2; %number of benchmark trials to compute average
  lbound = 0.01*(N*N); %default lowerbound on nnz
  bound = 10*10^6; %default upperbound on nnz
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
  ER = erdos_generator(N,prob); 
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
  for prob = p_min:grain:p_max
    ER = erdos_generator(N,prob,0,0);
    ER = ER(P,P);
    cut = getcut(ER,idxlist); 
    new_row = [nnz(ER) cut cut/nnz(ER) nnz(ER)/prod(size(ER))];
    lambdas = [lambdas; new_row];
  end
  disp(num2str(lambdas,'%.2f'))
end
