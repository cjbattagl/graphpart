%Input: Adjacency Matrix A, Weight parameter gamma
%Output: bisected adjacency matrix B = A(p,p)
function [B] = fennel (A, gamma, num_parts)
  vorder = randperm(size(A,1));

  m = size(A,1);
  e = nnz(A);

  alpha = e*((2^(gamma-1))/m^gamma);
  % experimental section alpha
%   alpha = sqrt(2)*e/m^1.5;

  parts = cell(num_parts);
  
  
  for v=vorder
    idxs = find(A(v,:));  %get N(v)

%     [~, loc] = ismember(idxs, p1); %set loc to idxs of intersect(N(v), S1)
%     loc = find(loc); %length(loc) is now number of edges between v, S1
%     p1newedges = length(loc)-length(idxs);  %p1newedges is number of edges between v, V \ S1 = e(v, V \ S1)
%     
    
    best_part_idx = 0;
    best_part_score = -inf;
    for p=1:num_parts
      [~, loc] = ismember(idxs, parts{p});
      loc = find(loc);
      p1newedges =  length(loc)-length(idxs);

      % if h(S2) + h(S1 U {v}) > h(S1) + h(S2 U {v})
      
      if (p1newedges - dc(length(parts{p}),alpha, gamma)) > best_part_score
        best_part_idx = p;
        best_part_score = p1newedges - dc(length(parts{p}),alpha, gamma);
      end
    end
    parts{best_part_idx} = [parts{best_part_idx} v];
  end
%   disp(size(p1))
%   disp(size(p2))
  
  B_ind = [];
  for i=parts
    B_ind = [B_ind  i{:}];
  end
  B = A(B_ind, B_ind);
%   B = A([p1 p2],[p1 p2]); 
end

function [y] = c(x, alpha, gamma)
  y = alpha * x^(gamma);
end

function [dy] = dc(x, alpha, gamma)
  dy = c(x+1, alpha, gamma) - c(x, alpha, gamma);
end
