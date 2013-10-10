% This is a benchmarking experiment where we grow ER(N,p) from t=0 -> infty
% We generate a directed graph with no self-edges
function [edges times zs] = grow_erdos(N,max_nnz,seed)
  
  iterations = 20; %number of iterations between p_min and p_max
  trials = 2; %number of benchmark trials to compute average
  lbound = 100; %default lowerbound on nnz
  bound = 1*10^6; %default upperbound on nnz
  times = []; %vector of benchmark results
  edges = []; %vector of nnz's
  zs = []; %average degrees
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

  for prob = p_min:grain:p_max
    ER = erdos_generator(N,prob);
    t = 0;

    for i=1:trials
      tic;
      ER*x;
      t = t + toc;
    end

    t = t / trials;
    times = [times t];
    edges = [edges nnz(ER)/(10^6)]; %megabytes of nnz
    zs = [zs prob*(N-1)];
  end

end
