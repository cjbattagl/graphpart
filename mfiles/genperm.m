function [P idxlist] = genperm(map)
  P = [];
  nparts = max(map(:,1));
  idxlist = zeros(nparts,2); %contains indices of partition matrix
  psizes = zeros(nparts,1);
  for i=1:nparts
    [t idx] = find(map(:,1)==i);
    verts = map(t,2);
    lo = length(P) + 1;
    P = [P; verts];
    hi = lo + length(t) - 1;
    idxlist(i,:) = [lo hi];
    psizes(i) = length(t);
  end
