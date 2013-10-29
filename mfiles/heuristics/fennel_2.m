%Input: Adjacency Matrix A, Weight parameter gamma
%Output: bisected adjacency matrix B = A(p,p)
function [B] = fennel_2 (A, gamma)
  vorder = randperm(size(A,1));

  m = size(A,1);
  e = nnz(A);

  alpha = e*((2^(gamma-1))/m^gamma);
  % experimental section alpha
%   alpha = sqrt(2)*e/m^1.5;

  p1 = [];
  p2 = [];
  
  for v=vorder
    idxs = find(A(v,:));  %get N(v)

    [~, loc] = ismember(idxs, p1); %set loc to idxs of intersect(N(v), S1)
    loc = find(loc); %length(loc) is now number of edges between v, S1
    p1newedges = -length(idxs) + length(loc); %p1newedges is number of edges between v, V \ S1 = e(v, V \ S1)
    clear loc
    
    [~, loc] = ismember(idxs, p2);
    loc = find(loc);
    p2newedges = -length(idxs) + length(loc);

    % if h(S2) + h(S1 U {v}) > h(S1) + h(S2 U {v})
    if (p1newedges - dc(length(p1),alpha, gamma)) >= (p2newedges - dc(length(p2),alpha,gamma))
      p1 = [p1 v];
    else
      p2 = [p2 v];
    end
  end
%   disp(size(p1))
%   disp(size(p2))
  
  
  B = A([p1 p2],[p1 p2]); 
end

function [y] = c(x, alpha, gamma)
  y = alpha * x^(gamma);
end

function [dy] = dc(x, alpha, gamma)
  dy = c(x+1, alpha, gamma) - c(x, alpha, gamma);
end
