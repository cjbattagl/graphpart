% Return A, adj. matrix for undirected E-R(N,p). If isDirected argument
% is set to true, return directed graph instead
% Existing ER generators for MATLAB aren't scalable for massive graphs because
% they create a dense random matrix.
% This implementation scales by initializing one row at a time.
function [A] = erdos_generator(N, p, isDirected, seed)
  if (nargin == 4)
    rand('seed',seed);
  end

  A = sparse(N,N);
  
  for i=1:N
    rand('seed',seed+i);
    vec = rand(N,1);
    vec = vec < p;
    vec = sparse(vec);
    A(:,i) = vec;
  end

  A = A(~eye(N)); %remove self-edges
  if (nargin == 3)
    if (~isDirected)
      T = tril(A);
      A = T + T';
    end
  else %assume undirected
    T = tril(A);
    A = T + T';
  end
end
