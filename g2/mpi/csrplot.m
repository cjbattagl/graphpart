  clf
  %A=load('mygraphtup.mat');
  A = load('mygraphpcsr.mat');
  %A = load('out_pcsr02.mat');

  %C=load('mygraphperm.mat');
  %D=load('mygraphpcsr.mat');
  A = spconvert(A);
  A(max(size(A)),max(size(A))) = 1;
  %B = spconvert(B);
  %C = spconvert(C);
  %D = spconvert(D);
  hold off;
  spy(A,'.');