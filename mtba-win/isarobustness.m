function rob = isarobustness(normed_data, row_scores, col_scores)
%Calculates the robustness of ISA biclusters
%
% Usage 
% >> rob = isarobustness(normed_data, row_scores, col_scores)
%
% Inputs
%   normed_data     -   the normalized input data, calculated with
%                       isaNormalize
%   row_scores      -   scores of row components of biclusterss
%   col_scores      -   scores of column components of biclusters
% 
% Output
%   rob             -   vector of robustness score of each bicluster
%
% See Also: isaFilterRobust
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if size(row_scores,2)~=size(col_scores,2)
  error('Row and column dimension mismatch');
end

Ec = normed_data.Ec;
Er = normed_data.Er;

for i=1:size(row_scores,2)
  row_scores(:,i) = row_scores(:,i)./sqrt(sumsqr(row_scores(:,i)));
end
for i=1:size(col_scores,2)
  col_scores(:,i) = col_scores(:,i)./sqrt(sumsqr(col_scores(:,i)));
end

if normed_data.hasNaN
  tempEr = Er;
  tempEr(isnan(tempEr)) = 0;
  tempEc = Ec;
  tempEc(isnan(tempEc)) = 0;
  rob1 = sum(bsxfun(@times, tempEr*row_scores, col_scores));
  rob2 = sum(bsxfun(@times, tempEc*col_scores, row_scores));
else
  rob1 = sum(bsxfun(@times, Er*row_scores, col_scores));
  rob2 = sum(bsxfun(@times, Ec*col_scores, row_scores));
end

rob1(rob1<0) = 0;
rob2(rob2<0) = 0;

rob = sqrt(rob1).*sqrt(rob2);
end