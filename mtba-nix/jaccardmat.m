function mat = jaccardmat(rows, cols)
%Jaccard similarity matrix function
%
% Usage 
% >> mat = jaccardmat(rows, cols)
%
% Inputs
%   rows  - matrix containing rows of biclusters
%   cols  - matrix containing columns of biclusters
%
% Output
%   Similarity matrix containg jaccard index between all biclusters
%
% See Also: jaccardind 
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

  le = size(rows,2);
  mat = zeros(le,le);

  for i=1:(le-1)
    for j=(i+1):le
      alle1 = rows(:,i) * cols(i,1)';
      alle2 = rows(:,j) * cols(j,1)';
      alle = alle1 + alle2;
      mat(i,j) = sum(sum(alle>1))/sum(sum(alle>0));
    end
  end
end
