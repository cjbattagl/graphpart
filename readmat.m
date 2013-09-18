% This function reads a `mat' style matrix file (graph)
% and returns it as a sparse matrix

function [B] = readmat(fname)
  B = load(fname);
  B = B.Problem.A;
end
