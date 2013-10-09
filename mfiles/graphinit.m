% This adds the paths to scripts I was using...
% CSparse is needed for cspy
% metismex is needed for metispart
% patohmat for patohmat ^_^

function graphinit

root = pwd;

cd external/matlab-bgl
path(path,pwd);
cd(root);
cd external/CXSparse/MATLAB/CSparse
disp('I hope you built CXSparse');
path(path,pwd);
cd(root);

%  cd CSparse/MATLAB/CSparse
%  path(path,pwd);
%  cd ~/metis-5.0.2/metismex
%  path(path,pwd);
%  cd ~/kron
%  cd ~/kron/patohmat
%  path(path,pwd);
%   cd ~/kron
end
