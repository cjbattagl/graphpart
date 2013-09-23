% creates a degree-likelihood plot of A
function [mean sd] = plotdegree(A)
  [a idx] = sort(sum(A));
  n = size(A,1);
  degs = full(sum(A));
  mean = sum(degs)/n;
  sd = std(degs);
  degs = degs(idx); %sort degrees
  bins = unique(degs);
  [bincounts] = histc(degs,bins);
  binprobs = bincounts ./ n;
  loglog(bins,binprobs);
end
