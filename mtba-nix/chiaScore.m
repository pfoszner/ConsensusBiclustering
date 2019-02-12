function [Tscore1 Bscore1 Tscore2 Bscore2 SBscore TSscore] = ...
                                        chiaScore(x, biClustResult, number)
%Compute Chia and Karuturi scores
% BK Chia & RK Karuturi
% Differential co-expression framework to quantify goodness of biclusters 
% and compare biclustering algorithms
% Algorithms for Molecular Biology, 5, 23. 
%
% Usage
% >> [Tscore1 Bscore1 Tscore2 Bscore2 SBscore TSscore] = ...
%                                       chiaScore(x, biClustResult, number)
%
% *Inputs*
%   x             - data matrix
%   biClustResult - bicluster result
%   number        - cluster number/index for computing scores
%
% *Outputs*
%   Tscore1   - Row effects for chosen bicluster(within)
%   Bscore1   - Column effect for chosen bicluster(within)
%   Tscore2   - Row effect for chosen bicluster(outside)
%   Bscore2   - Column effect for chosen bicluster(outside)
%   SBscore   - Ranking Score
%   TSscore   - Stratification Score
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

  biclCols = biClustResult.Clust(number).cols;
  biclRows = biClustResult.Clust(number).rows;
  xmat = x(biclRows,biclCols);

  yBar_j = mean(xmat,1); % Column Mean
  yBar_i = mean(xmat,2); % Row Mean
  x_vec  = xmat(:);
  yBar   = mean(x_vec);
  E = 0;

  for i=1:length(biclRows)
    for j=1:length(biclCols)
      E = E + (xmat(i,j) - yBar_i(i) - yBar_j(j) + yBar)^2;
    end
  end
  E = E/((size(xmat,2) - 1) * (size(xmat,1) - 1));

  T1 = sum(sum((yBar_i).^2))/size(xmat,1) - E/size(xmat,2);
  B1 = sum(sum((yBar_j).^2))/size(xmat,2) - E/size(xmat,1);

  xmat2 = x(biclRows, :);
  xmat2(:,biclCols) = [];

  yBar_j2 = mean(xmat2, 1); % Column Mean
  yBar_i2 = mean(xmat2, 2); % Row Mean
  x_vec2  = xmat2(:);
  yBar    = mean(x_vec2);
  E2 = 0;

  for i=1:size(xmat2,1)
    for j=1:size(xmat2,2)
      E2 = E2 + (xmat2(i,j) - yBar_i2(i) - yBar_j2(j) +ybar2)^2;
    end
  end
  E2 = E2/((size(xmat2,2) - 1) * (size(xmat2,1) - 1));

  T2 = sum(sum(yBar_i2.^2))/size(xmat2,1) - E/size(xmat2,2);
  B2 = sum(sum(yBar_j2.^2))/size(xmat2,2) - E/size(xmat2,1);

  a = 0.01;

  SB = log(max([T1+a B1+a]) / max([T2+a B2+a]));
  if SB>0, TS = log((T1+a)/(B1+a)); else TS = log((T2+a)/(B2+a)); end

  Tscore1 = T1, Bscore1 = B1, Tscore2 = T2, Bscore2 = B2;
  SBscore = SB, TSscore = TS;
end
