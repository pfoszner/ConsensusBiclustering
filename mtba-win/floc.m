function biClustResult = floc(data, numBiclust, pRow, pColumn, resTh, ...
                                minRow, minCol, nIter, blocRow, blocColumn)
%FLexible Overlapped biClustering Algorithm
% Yang J, et al. 
% An improved biclustering method for analyzing gene expression profiles. 
% Int. J. Artif. Intell. T. 2005;14:771-790.
%
% Inputs
%   data        -   input data matrix
%   numBiclust  -   the number of biclusters searched
%   pRow        -   genes initial probability of membership to the 
%                   biclusters
%   pColumn     -   samples initial probability of membership to the 
%                   biclusters
%   resTh       -   residue threshold
%   minRow      -   minimal number of gene per bicluster
%   minCol      -   minimal number of conditions per bicluster
%   nIter       -   number of iterations
%   blocRow     -   a matrix indicating the directed initialisation for
%                   the genes
%   blocColumn  -   a matrix indicating the directed initialisation for
%                   the conditions
%
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] if Row i is in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] if Col j is in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices.
%       MatRes      - A matrix describing the biclusters with the following
%                     columns:
%                     Residue, Volume, Genes(rows), Conditions(columns), 
%                     Row Variance
%
% See Also: RESIDU
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<2
  numBiclust = 20;
end        
if nargin<3
  pRow = 0.5;
  pColumn = pRow;
end
if nargin<4
  pColumn = pRow;
end
if nargin<5
  resTh = [];
end
if nargin<6
  minRow = 8;
end
if nargin<7
  minCol = 6;
end
if nargin<8
  nIter = 500;
end
if nargin<9
  blocRow = [];
end
if nargin<10
  blocColumn = [];
end


if isempty(resTh), resTh = residu(data)/10; end

vecData = data(:);
nbRows = size(data,1);
nbCols = size(data,2);

if ~isempty(blocRow) || ~isempty(blocColumn)
  numBiclust = max([size(blocRow,1) size(blocColumn,1) numBiclust]);
end

vecBicRow = zeros(nbRows, numBiclust);
vecBicCol = zeros(nbCols, numBiclust);

if isempty(blocRow), blocRow = zeros(nbRows, numBiclust);
else vecBicRow(:,1:size(blocRow,1)) = blocRow; end
vecBlocRow = vecBicCol(:);
vecBicRow = vecBicRow(:);

if isempty(blocColumn), blocColumn = zeros(nbCols, numBiclust);
else vecBicCol(:,1:size(blocColumn,1)) = blocColumn; end
vecBlocCol = vecBicCol(:);
vecBicCol = vecBicCol(:);

% Generate Random
rand1 = rand(1, numBiclust*nbRows);
rand2 = rand(1, numBiclust*nbCols);

vecBicRow(rand1 < pRow) = 1;
vecBicCol(rand2 < pColumn) = 1;

vecResvol_bic = zeros(numBiclust*4,1);

mfloc( vecData, nbRows, nbCols, vecBicRow, vecBicCol, ...
                  vecResvol_bic, resTh, numBiclust, minRow, minCol,...
                  nIter, vecBlocRow, vecBlocCol);
bicRow = reshape(vecBicRow, nbRows, numBiclust)'; %Row wise reshape
bicCol = reshape(vecBicCol, nbCols, numBiclust)'; %Row wise reshape

matResvol_bic = zeros(numBiclust, 5);
tempmat = reshape(vecResvol_bic, 4, numBiclust)'; %Row wise reshape
matResvol_bic(:,1:4) = tempmat(:,1:4);
for i=1:numBiclust
  v = var(data(find(bicRow(i,:)==1), find(bicCol(i,:)==1)),0,2);
  matResvol_bic(i,5) = mean(v(:));
end
% Columns of matResvol_bic
% --------------------------------------------------------
% | Residue | Volume | Genes | Conditions | Row Variance |
% --------------------------------------------------------

% Result
RowxNum = bicRow';
NumxCol = bicCol;
for i=1:numBiclust
 rows = find(RowxNum(:,i)>0);
 cols =find(NumxCol(i,:)>0);
 cluster(i) = struct('rows', rows, 'cols', cols);
end

% Returning result
biClustResult = struct('RowxNum',RowxNum, 'NumxCol',NumxCol, ...
                        'ClusterNo', numBiclust, 'Clust', cluster, 'MatRes', matResvol_bic);
end

