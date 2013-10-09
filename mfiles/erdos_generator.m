% Return A, adj. matrix for undirected E-R(N,p). If isDirected argument
% is set to true, return directed graph instead
function [A] = erdos_generator(N, p, isDirected)
  A = (rand(N) <= p);
  A = A(~eye(N)); %remove self-edges
  if (nargin == 3)
    if (~isDirected)
      T = tril(A);
      A = T + T';
    end
  else %assume directed
    T = tril(A);
    A = T + T';
  end
end
