function plotclust(res, x, legende, noC, Title)
%Plot comparing values inside biclusters to values outside
%
% Inputs
%   res       - biClustResult
%   x         - data matrix for which bicluster was computed
%   legende   - `true` to display legend.
%   noC       - Number of Clusters to plot. Default:5
%   Title     - String to be displayed as title.
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<3
  legende = false;
  noC = 5;
  Title = 'Plotclust';
end

for i=1:min([res.ClusterNo, noC])
  identq = res.NumxCol(i,:);
  identper = res.RowxNum(:,i);
  
  a = sum(identper(:));
  le = sum(identq(:));
  
  if a==0, break; end
  
  
  mat = zeros(le, 3);
  mat(:,1) = mean(x(identper,identq),1);
  mat(:,2) = mean(x(~identper,identq),1);
  mat(:,3) = median(x(~identper,identq),1);
  
  subplot(noC,1,i); % Look for subplotting
  bar(mat); title({Title; ['Cluster:' num2str(i) ', Size:' num2str(a)]});
  if legende
    legend('Value', 'Mean', 'Median')  ;
  end
end
end
