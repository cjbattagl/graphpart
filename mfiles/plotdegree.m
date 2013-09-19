% creates a degree-likelihood plot of A
function [] = plotdegree(A)
  [a idx] = sort(sum(A));
  n = size(A,1);
  degs = full(sum(A));
  degs = degs(idx); %sort degrees
  bins = unique(degs);
  [bincounts] = histc(degs,bins);
  binprobs = bincounts ./ n;
  loglog(bins,binprobs);
end
