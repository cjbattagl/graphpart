% This was a function I threw together to implement FENNEL
% It's gross and hard to read right now and I don't even remember if it works
% I will fix it X_X

function [B] = fennel (A, gamma)
  vorder = bfsorder(A);

  m = size(A,1);
  e = nnz(A);

  alpha = e*((2^(gamma-1))/m^gamma);

  p1 = [];
  p2 = [];
  gp1 = 0; % maintain e(S1,V\S1)
  gp2 = 0; % maintain e(S2,V\S2)

  for v=vorder
    idxs = find(A(v,:));  %get N(v)
    [tf, loc] = ismember(idxs, p1); %set loc to idxs of intersect(N(v), S1)
    loc = find(loc); %length(loc) is now number of edges between v, S1
    p1newedges = length(idxs) - length(loc); %p1newedges is number of edges between v, V \ S1 = e(v, V \ S1)
    tempgp1 = gp1 + p1newedges; %tempgp1 is e(S1 U {v}, V \ (S1 U {v}))
    g1 = tempgp1 - c(length(p1), alpha, gamma); %h(S1 U {v})
    
    [tf, loc] = ismember(idxs, p2);
    loc = find(loc);
    p2newedges = length(idxs) - length(loc);
    tempgp2 = gp2 + p2newedges;
    g2 = tempgp2 - c(length(p2), alpha, gamma); %h(S2 U {v})

    % if h(S2) + h(S1 U {v}) > h(S1) + h(S2 U {v})
    if (((gp2 - c(length(p2),alpha,gamma))+ g1) > ((gp1- c(length(p1),alpha,gamma)) + g2))
      p1 = [p1 v];
      gp1 = g1;
    else
      p2 = [p2 v];
      gp2 = g2;
    end
  end
  %disp(size(p1))
  %disp(size(p2))
  B = A([p1 p2], [p1 p2]); 
end

function [y] = c(x, alpha, gamma)
  y = alpha * x^gamma;
end

