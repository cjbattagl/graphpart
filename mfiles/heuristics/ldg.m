function [B] = ldg(A)
  N = size(A,1);
  vorder = randperm(N);
  e = nnz(A);
  eps = N/(10);
  C = (N+eps)/2;

  p1 = [vorder(1)];
  p2 = [];

  for v=vorder(2:N)
    idxs = find(A(v,:));  %get N(v)
    [tf, loc] = ismember(idxs, p1); %set loc to idxs of intersect(N(v), S1)
    loc = find(loc); %length(loc) is now number of edges between v, S1
    e_v_s1 = length(loc) + 1;
    
    [tf, loc] = ismember(idxs, p2);
    loc = find(loc);
    e_v_s2 = length(loc) + 1;
    
    w_t_1 = 1-((length(p1))/C);
    w_t_2 = 1-((length(p2))/C);

    if ((w_t_1 * e_v_s1) > (w_t_2 * e_v_s2))
      p1 = [p1 v];
    else
      p2 = [p2 v];
    end
  end
  disp(size(p1))
  disp(size(p2))
  B = A([p1 p2], [p1 p2]); 
end
