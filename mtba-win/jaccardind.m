function index = jaccardind(bicres1, bicres2)
%Adaption of Jaccard Index for clustering
%
% Usage
% >> index = jaccardind(bicres1, bicres2)
%
% Inputs
%   bicres1, bicres2 - bicluster result structure
% 
% Output
%   Jaccard index value for the clusters
%
% See Also: jaccardmat
% 
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India
 
  index = jaccard1(bicres1, bicres2)/max([jaccard1(bicres1, bicres1) jaccard1(bicres2, bicres2)]);
end

function res = jaccard1(bicres1, bicres2)
  le1 = bicres1.ClusterNo;
  le2 = bicres2.ClusterNo;
  jacvec = [];

  for i=1:le1
    jacvec2 = 0;
    for j=1:le2
      alle1 = bicres1.RowxNum(:,i) * bicres1.NumxCol(i,:)';
      alle2 = bicres2.RowxNum(:,j) * bicres2.NumxCol(j,:)';
      alle = alle1 + alle2;
      loalle = alle>0;
      loalle1 = alle1>0;
      loalle2 = alle2>0;
      jacvec2 = jacvec2 + (sum(loalle1(:)) + sum(loalle2(:)) - sum(loalle(:)))...
                          /sum(loalle(:));
    end
    jacvec(end+1) = jacvec2;
  end

  res = sum(jacvec(:))/max([le1 le2]);
end


         
         
