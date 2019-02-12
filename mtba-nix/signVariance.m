function v = signVariance(x, resultSet, number, dimension)
% Preliminary measure of sign coherence of a bicluster
%
% S. C. Madeira, A. L. Oliveira
% Biclustering algorithms for biological data analysis: a survey.
% IEEE/ACM Trans Comput Biol Bioinform. 2004 Jan-Mar ;1(1):24-45.
%
% Inputs
%   x         - The data matrix from which biclusters were identified
%   resultSet - The bicluster result calculated from the data matrix
%   number    - The number(index) of bicluster within the result.
%   dimension - 'both': overall variance
%               'row' : row variance
%               'col' : column variance
% Output
%   v   - Corresponding variance of genes/conditions as the average of the
%   sum of euclidean distances between all rows and/or columns of the
%   bicluster.
%
% Note: The lower the value returned, the more coherent the bicluster is.
%
% See Also: variance
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<4
  dimension = 'both';
end

rows = resultSet.Clust{number}.rows;
cols = resultSet.Clust{number}.cols;

A = x(rows,cols);

if strcmp(dimension, both)
  rv = variance(sig(A,2));
  cv = variance(sig(A),2);
  [n m] = size(A);
  
  v = (n.*rv+m.*cv)./(n+m);
else
  if strcmp(dimension, 'row'), v = variance(sig(A,2),2);  % Same as Matlab way
  elseif strcmp(dimension, 'col'), v = variance(sig(A),1);
  end
end
end

function ret = sig(x ,rowOrColumn)
% Similar to div() & dif() but only change of slope calculated
%
% *Inputs*
%   x           - matrix (n x m)
%   rowOrColumn - 1: sign by columns
%                 2: sign by rows
%                 Default:1
%
% *Output*
%   matrix ((n-1) x m) or (n x (m-1)) (depending on rowOrColumn)
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
    for i=1:n
      if x(i,j)>x(i,j-1), ret(i,j)=1; end
      if x(i,j)<x(i,j-1), ret(i,j)=-1; end
      if x(i,j)==x(i,j-1), ret(i,j)=0; end
    end
  end
elseif rowOrColumn==2
  ret = NaN(n,m);
  ret(1,:) = 0;
  for i=2:n
    for j=1:m
      if x(i,j)>x(i-1,j), ret(i,j)=1; end
      if x(i,j)<x(i-1,j), ret(i,j)=-1; end
      if x(i,j)==x(i-1,j), ret(i,j)=0; end
    end
  end
end
end