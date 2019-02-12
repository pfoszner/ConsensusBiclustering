function v = additiveVariance(x, resultSet, number, dimension)
% Preliminary measure of additive coherence of a bicluster
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
%
% Output
%   v   - Corresponding variance of genes/conditions as the average of the
%   sum of euclidean distances between all rows and/or columns of the
%   bicluster.
%
% See Also: dif
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<4
  dimension = 'both';
end

rows = resultSet.Clust(number).rows;
cols = resultSet.Clust(number).cols;

A = x(rows,cols);

if strcmp(dimension, both)
  rv = variance(dif(A,2));
  cv = variance(dif(A),2);
  [n m] = size(A);
  
  v = (n.*rv+m.*cv)./(n+m);
else
  if strcmp(dimension, 'row'), v = variance(dif(A,2),2);  % Same as Matlab way
  elseif strcmp(dimension, 'col'), v = variance(dif(A),1);
  end
end
end

