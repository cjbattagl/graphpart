% Return the permutation of A that results from a BFS starting from node i
% (Following a directed edge from vertex u to vertex v if A(u,v) != 0)
% Note that not all nodes may be reachable

function [B explored] = bfsorder(A, i)
explored = [i];
frontier = [];
terminate = 0;

while (~terminate)
   [row col] = find(A(explored,:));
   if (size(col,1) > 1) col = col'; end;
   frontier = uniqueorder(col);
   num_explored = length(explored);
   explored = uniqueorder([explored frontier]);
   if (length(explored) == num_explored)
     terminate = 1;
   end
end
B = A(explored,explored);
end

%remove duplicates without changing the order of X
function [Y] = uniqueorder(X)
  [Xs, SortVec] = sort(X);
  UV(SortVec) = ([1 diff(Xs)] ~= 0);
  Y = X(UV);
end