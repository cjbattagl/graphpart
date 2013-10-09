% This adds the paths to scripts I was using...
% CSparse is needed for cspy
% metismex is needed for metispart
% patohmat for patohmat ^_^

function graphinit

curr = pwd;
cd external/mit_matlab_tools
path(path,pwd);
cd curr;

%  cd CSparse/MATLAB/CSparse
%  path(path,pwd);
%  cd ~/metis-5.0.2/metismex
%  path(path,pwd);
%  cd ~/kron
%  cd ~/kron/patohmat
%  path(path,pwd);
%   cd ~/kron
end
