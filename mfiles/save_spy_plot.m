% This function saves a spy plot of the inputted matrix as a png
% with file name fname, in the current working directory

function save_spy_plot(A, fname)
  disp('Saving...');
  [s,M,H] = cspy(A,1024);
  imwrite(M,H,fname);
end
