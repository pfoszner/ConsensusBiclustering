function ret = dif(x, rowOrColumn)
%Equivalent to div(), uses subtraction rather than division between
%each row/column and the previous row/column
%
% Inputs
%   x           - matrix (n x m)
%   rowOrColumn - 1: subtraction by columns
%                 2: subtraction by rows
%                 Default:1
%
% Output
%   matrix ((n-1) x m) or (n x (m-1)) (depending on rowOrColumn)
%
% See Also: additiveVariance 
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<2
  rowOrColumn = 1;
end

[n m] = size(x);

if rowOrColumn==1 % Column
  ret = NaN(n,m);
  ret(:,1) = 0;
  for j=2:m
    ret(:,j) = x(:,j) - x(:,j-1);
  end
elseif rowOrColumn==2
  ret = NaN(n,m);
  ret(1,:) = 0;
  for i=2:n
    ret(i,:) = x(i,:) - x(i-1,:);
  end
end
end