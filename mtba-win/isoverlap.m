function olap = isoverlap(biClustResult)
% Check if biclusters are overlapping
%
% See Also: filterclust
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

rows = max(max(sum(biClustResult.RowxNum,2))); % Row Sum
cols = max(max(sum(biclustResult.NumxCol,1))); % Column Sum

if rows>1 && cols>1
  disp('There are overlapping rows and columns in the Bicluster result.');
  res = true;
elseif rows>1
  disp('There are overlapping rows in the Bicluster result');
  res = true;
elseif cols>1
  disp('There are overlapping columns in the Bicluster result');
  res = true;
else
  disp('There are no overlapping rows or columns in the Bicluster result');
  res = false;
end

olap = struct('Overlapping',res,'MaxbiclusterRows',rows,...
              'MaxbiclusterCols',cols);
  
end