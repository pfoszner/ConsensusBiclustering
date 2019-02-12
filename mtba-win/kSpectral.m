function biClustResult = kSpectral(matrix, normalization, numberEigenvalues, minr, minc, withinVar)
%Spectral Biclustering 
% Kluger Y., Basri R., Chang J.T. & Gerstein M., 
% Spectral Biclustering of Microarray Data: Coclustering Genes and Conditions
% Genome Research 2003. 
%
% Usage 
% >> biClustResult = kSpectral(matrix, normalization, numberEigenvalues,...
%                               minr, minc, withinVar)
%
% Inputs:
%   matrix            - data matrix
%   normalization     - Normalization method to be applied to the matrix. As
%                       explained in Kluger et al, three methods are allowed:
%                       1)"log":Logarithmic Transformation
%                       2)"irrc":Independent Rescaling of Rows and Columns
%                       3)"bistochast":Bistochastization
%                     If "log" normalization is used, be sure you can apply
%                     logarithm to elements in data matrix, if there are 
%                     values under 1, it automatically will sum to each 
%                     element in mat (1+abs(min(mat))) Default is "log", 
%                     as recommended by Kluger et al.
%   numberEigenvalues - number of eigenvalues considered to find biclusters.
%                      Each row eigen-vector combines with all column 
%                      eigen-vectors for the first 'numberEigenvalues'
%                      eigenvalues.
%                      Note: a high number could decrease dramatically time
%                      performance. Usually, only the very ﬁrst 
%                      eigenvectors are used. With "irrc" and 
%                      "bistochastization" methods, ﬁrst eigenvalue 
%                      contains background (irrelevant) information, so it 
%                      is ignored.
%   minr              - Minimum number of rows a bicluster must have.
%   minc              - Minimum number of columns a bicluster must have.
%   withinVar         - maximum within variation allowed
% 
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] if Row i is in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] if Col j is in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices.
%
% See Also: normalize
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

%##########################################################################
%       Error Checking and Default Values
%##########################################################################
if nargin<1
  error('input :  No matrix as input');
end

if nargin<2
  normalization='log';
  numberEigenvalues=3;
  minr=2;
  minc=2;
  withinVar=1;
end

A = matrix;
n = size(A,1);
m = size(A,2);
discardFirstEigenvalue = true;

% 1) Normalization
if strcmp(normalization,'log')
  if (min(A(:))<1)
    A = bsxfun(@plus,A,1+abs(min(A(:))));
  end
  K = normalize(A,'logt') ; 
  discardFirstEigenvalue = false;
else
  if numberEigenvalues==1
    warning('Only 1 eigenvalue selected with a normalization that gives first eigenvalue as background. Increasing to 2 eigenvalues');
    numberEigenvalues = 2;
  end
  
  Ab = bsxfun(@plus,A,abs(min(A(:))));
  if strcmp(normalization,'irrc')
    K = normalize(Ab,'irrc'); 
  elseif strcmp(normalization,'bistochast')
    K = normalize(Ab,'bistochast'); 
  end
end

% 2) SVD Decomposition
[U,~,V] = svd(K,'econ');

% 3) Vector processing, basically k-means
result = postprocess(U,V,numberEigenvalues,2,min([n/minr,100]),2,m/minc);

% 4) Taking all possible eigenvector combination biclusters
srows = {};
scols = {};
init=1;
if(discardFirstEigenvalue)
  init=2;
end

for k=init:numberEigenvalues
  init2=k;
  for l=init2:numberEigenvalues
    numGenes = result.numgenes(k);
    numExpr  = result.numexp(l);
    
    for i=1:numGenes
      for j=1:numExpr
        % Some hocus pocus, will get error here for sure
        indexArow=1:size(A,1);
        indexAcol=1:size(A,2);
        srows{end+1} = indexArow(result.eigenGeneCluster(k,:)==i);
        scols{end+1} = indexAcol(result.eigenexprCluster(l,:)==j);
      end
    end
  end
end

% 5) Discard non relevant biclusters
srowsOK = {};
scolsOK = {};
ss = [];
%   ret = [];

for i=1:size(srows,2)
  ncols = size(scols{i},2);
  nrows = size(srows{i},2);
  wv    = withinvar(A(srows{i},scols{i}), ncols, nrows);
  
  if (wv<withinVar && ncols>minc && nrows>minr)
    srowsOK{end+1} = srows{i};
    scolsOK{end+1} = scols{i};
    ss = [ss, wv];
  end
end

if size(srowsOK)==0
  warning('No biclusters found');
  biClustResult = struct('RowxNum',NaN, 'NumxCol',NaN, ...
    'ClusterNo', 0, 'Clust', NaN);
else
  RowxNumber = false(size(A,1), size(srowsOK,2));
  ColxNumber = false(size(A,2), size(srowsOK,2));
  
  for i=1:size(srowsOK, 2)
    RowxNumber(srowsOK{i},i) = true;
  end
  for i=1:size(scolsOK, 2)
    ColxNumber(scolsOK{i},i) = true;
  end
  for i=1:size(srowsOK, 2)
    rows = find(RowxNumber(:,i)>0);
    cols =find(ColxNumber(:,i)>0);
    cluster(i) = struct('rows', rows, 'cols', cols);
  end
  % Returning result
  biClustResult = struct('RowxNum',RowxNumber, 'NumxCol',ColxNumber', ...
    'ClusterNo', size(srowsOK,2), 'Clust', cluster);
end
end

%
%  Within variation of a matrix by rows
%  Computes the row mean and then the euclidean distance of each row to the mean.

function within = withinvar(x,n,m)
if n==1 % One row
  within = 0;
  
elseif m==1 % One column
  centroid = mean(x(:));
  distances = sqrt(sumsqr(bsxfun(@minus,x,centroid)));
  within = sum(distances)/n;
else
  centroid = mean(x,1);
  distances = sqrt(sum(((bsxfun(@plus,-x,centroid))).^2,1));
  within = (sum(distances))/n;
end
end
