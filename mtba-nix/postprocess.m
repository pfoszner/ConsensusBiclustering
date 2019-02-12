function dev = postprocess(u,v,maxeigen,minCG,maxCG,minCE,maxCE)
%Postprocessing of eigenvectors as described in:
% Kluger, Y.; Basri, R.; Chang, J.T. & Gerstein, M.,
% Spectral Biclustering of Microarray Data: Coclustering Genes and Conditions
% Genome Research 2003.
%
% Apply iterative kmeans clustering to all possibel combinations of
% eigenvectors  of SVD analysis

% Usage
% >> dev = postprocess(u,v,maxeigen,minCG,maxCG,minCE,maxCE)
%
% Inputs:
%   u,v         - decomposition result from svd
%   maxeigen    - maximum number of eigenvectors of dec processed. 
%   minCG       - minimum number of clusters in which eigengenes will be clustered.
%   maxCG       - maximum number of clusters in which eigengenes will be clustered
%   minCE       - minimum number of clusters in which eigenexpressions will be clustered. 
%   maxCE       - maximum number of clusters in which eigenexpressions will be clustered
%
% Outputs:
%
% See Also: kSpectral
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India


n=size(u,1);
m=size(v,1);

max1=min([n m]);
if maxeigen>max1
  warning('Number of required eigenvectors exceeds freedom degrees');
  maxeigen = max1;
end

dev = struct('eigenGeneCluster',NaN(maxeigen,n),'eigenexprCluster',...
  NaN(maxeigen,m),'numgenes',1:maxeigen,'numexp',1:maxeigen);

for i=1:maxeigen
  for j=1:maxeigen
    u1 = u(:,i);
    v1 = v(:,j);
    % Best cluster y modified K-means
    clust = iterativeKmeans(u1,minCG,maxCG);
    dev.eigenGeneCluster(j,:) = clust;
    dev.numgenes(j) = max(clust(:));
    
    % Same for eigenarrays
    clust = iterativeKmeans(v1,minCE,maxCE);
    dev.eigenexprCluster(j,:) = clust;
    dev.numexp(j) = max(clust(:));
  end
end
end

function clust = iterativeKmeans(x, minimum, maximum)
% *Usage*  clust = iterativeKmeans(x, minimum, maximum)
%
% *Inputs*:
%   x       - input matrix
%   minimum - minimum number of clusters
%   maximum - maximum number of clusters
%
% *Outputs*:
%   clust   - vector  containing the cluster indices of each point

n = size(x,1);

if minimum<2
  error('Clustering must divide data in atleast 2 clusters. minimum value should be atleast 2');
end
if maximum>=n
  %error('Clustering must divide data in as much as n-1 clusters. maximum value should be less than size(x)-1');
end
if maximum<minimum
  error('maximum value should be greater than minimum');
end

numClusterings = maximum-minimum+1;
k = minimum;
ss = [];
while k<=maximum
  [~,~,sumd] = kmeans(x,k,'Options',statset('MaxIter',100),'emptyaction','singleton');
  ss = [ss; mean(sumd)];
  k = k+1;
end

difference = ss(2:size(ss))-ss(1:size(ss)-1);
md = mean(difference);
optim = 1;

for i=1:size(difference)
  if difference(i)>md
    optim = i;
    break;
  end
end
optim = optim+minimum-1;

clust = kmeans(x,optim,'Options',statset('MaxIter',100),'emptyaction','singleton');
end
