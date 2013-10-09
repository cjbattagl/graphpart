% Return A, adj. matrix for undirected E-R(N,p). If isDirected argument
% is set to true, return directed graph instead
% Existing ER generators for MATLAB aren't scalable for massive graphs because
% they create a dense random matrix.
function [A] = erdos_generator(N, p, isDirected, seed)
  if (nargin == 4)
    rand('seed',seed);
  end

  A = spones(sprand(N,N,p)); 
  
  if (nargin == 3)
    if (~isDirected)
      T = tril(A);
      A = T + T';
    end
  else %assume undirected
    T = tril(A);
    A = T + T';
  end

  %remove diagonals
  [nRows,nCols] = size(A);
  A(1:(nRows+1):nRows*nCols) = 0;
end
