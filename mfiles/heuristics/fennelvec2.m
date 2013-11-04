%This optimizes by recognizing that we only need to look in iteratively larger
%submatrices of A(vorder,vorder)
function [B] = fennelvec2 (A, gamma, num_parts)
  vorder = randperm(size(A,1));

  m = size(A,1);
  e = nnz(A);

  %alpha = e*((2^(gamma-1))/m^gamma);
  alpha = sqrt(2)*e/m^1.5;
  parts = zeros(num_parts,m,'logical');
  C = A(vorder,vorder);
  
  for v=1:m
    idxs = find(C(v,1:v));  %get N(v)    
    best_part_idx = 0;
    best_part_score = -inf;
    for p=1:num_parts
      p1newedges = nnz(parts(p,:)(idxs));

      if (p1newedges - dc(nnz(parts(p,:)),alpha, gamma)) > best_part_score
        best_part_idx = p;
        best_part_score = p1newedges - dc(nnz(parts(p,:)),alpha, gamma);
      end
    end
    parts(best_part_idx,v) = 1;
  end

  B_ind = [];
  for i=1:num_parts
    ind = find(parts(i,:));
    B_ind = [B_ind  ind];
  end

  B = C(B_ind, B_ind);
end

function [y] = c(x, alpha, gamma)
  y = alpha * x^(gamma);
end

function [dy] = dc(x, alpha, gamma)
  dy = c(x+1, alpha, gamma) - c(x, alpha, gamma);
end