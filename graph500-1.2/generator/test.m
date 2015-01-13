  A=load('mygraph.mat');
  A(:,1:2) = A(:,1:2) + 1;
  A = spconvert(A);
  spy(A,'.');
