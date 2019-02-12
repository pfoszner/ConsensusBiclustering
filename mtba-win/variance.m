function v = variance(x, rowOrColumn)
% Computes variance of input matrix x by rows or columns as the
% sum of euclidean distances between all rows/columns, divided by 1/n(n-1),
% being n the number of rows or columns.
% Zero variance implies an homogeneous matrix by rows or columns
%
% Note: This mesure is only applicable for biclusters based on constant
% values (by row, column or both) or coherent evolutions based on grouping
% next values 
%
% See Also: additiveVariance, overallVariance, constantVariance,
%           signVariance, multiplicativeVariance
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<2
  rowOrColumn = 1; %Column
end

[n m] = size(x);
v = 0;

if rowOrColumn==1 % Column
%   distan = squareform(pdist(x'));
  v = (1/(m*(m-1))) * sum(sum(pdist(x')));
elseif rowOrColumn==2 % Row
%   distan = squareform(pdist(x));
  v = (1/(n*(n-1))) * sum(sum(pdist(x)));
end
end