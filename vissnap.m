% Show a spy plot of the snap-generated graph with filename fname

function vissnap(fname)
  A = readsnap(fname);
  [a idx] = sort(sum(A));
  idx = fliplr(idx);
  A=A(idx,idx);
  cspy(A);
end