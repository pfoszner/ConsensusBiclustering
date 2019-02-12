function bicsout = filterclust(bicsin, overlapThreshold, maxNumber)
%list of the largest "<=maxNumber" biclusters with less overlap than "overlapThreshold"
%
% Usage
% >> bicsout = filterclust(bicsin, overlapThreshold, maxNumber)
%
% Input
%   bicsin            - cell array of biclusters (One way is the 
%                       biClusterResult.Clust)
%   overlapThreshold  - Percentage of overlap. Must be between [0,1]
%   maxNumber         - Default:100
%
% See Also: isoverlap
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<3
  maxNumber = 100;
end
if nargin<2
  overlapThreshold = 0.25;
end

cont = 0;

for i=1:length(bicsin)
  sizes(i) = length(bicsin(i).rows) * length(bicsin(i).cols);
end

bicsout = [];
[~,ord] = sortrows(sizes(:));
for i=bicsin{ord}
  if cont>=maxNumber, break; end
  insert = true;
  for j=bicsout
    if overlap(i,j)>=overlapThreshold
      insert = false;
      break;
    end
  end
    if insert
      bicsout{end+1}=i;
      cont = cont + 1;
    end
end
end

function ol = overlap(bic1, bic2)
% the amount of overlap between two biclusters, measured in terms of 
% percentage respect to bic1 (bic1 and bic2 are bicluster objects with 
% fields rows and cols)

  ol = length(intersect(bic1.rows, bic2.rows))/length(bic1.rows) * ...
       length(intersect(bic1.cols, bic2.cols))/length(bic1.cols);
end