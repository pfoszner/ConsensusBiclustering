function biClustResult = itersa(data,thr_row,thr_col,nSeeds, direction)
%Iterative signature algorithm by
% Sven Bergmann, Jan Ihmels and Naama Barkai
% Iterative signature algorithm for the analysis of large-scale gene expression data
% Phys Rev E Nonlin Soft Matter Phys
%
% Usage 
% >> biClustResult = itersa(data,thr_row,thr_col,nSeeds, direction)
%
% Inputs
%   data        - the input numeric matrix
%   thr_row     - a vector of row threshold parameters for which ISA runs
%   thr_col     - a vector of column thteshold parameters for which ISA
%                 runs
%   nSeeds      - Number of seeds to use. (Default:100)
%   direction   - String cell array of length two, one for rows, second for
%                 columns. Specifies whether we are interested in 
%                 rows/columns that higher('up') than average, lower 
%                 ('down') than average, or both ('updown').
%                 Default: {'updown','updown'}
%
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] 
%                     if Row i is in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] 
%                     if Col j is in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices_
%
% See Also: generateSeeds, isaFilterRobust, isaNormalize, isaUnique
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<1
  error('input :  No matrix as input');
end

if ~ismatrix(data)
  error('input : data must be a matrix');
end

if nargin<2
  thr_row   = 1:0.5:3;
  thr_col   = 1:0.5:3;
  nSeeds    = 100;
  direction = {'updown','updown'};
end

% Normalize the matrix
normed_data = isaNormalize(data);

% Generate seeds
row_seeds = generateSeeds(size(data,1), nSeeds);

% Determine thresholds
[tempx tempy] = meshgrid(thr_row, thr_col);
thr_list = [tempy(:) tempx(:)];

% Apply ISA, for all thresholds
for i=1:length(thr_list)
  isaresults{i} = isaIterate(normed_data, row_seeds, [], thr_list(i,1), thr_list(i,2), direction);
end

% Make it unique for every threshold combination
for i=1:length(isaresults)
  isaresults{i} = isaUnique(normed_data, isaresults{i});
end

% Filter according to robustness
for i=1:length(isaresults)
  isaresults{i} = isaFilterRobust(data, normed_data, isaresults{i}, 1, row_seeds, []);
end

% Merge them all
rows=[];
cols=[];
seeddata=[];
N=0;
for i=1:length(isaresults)
  rows = horzcat(rows, isaresults{i}.Rows);
  cols = horzcat(cols,isaresults{i}.Columns);
  seeddata = vertcat(seeddata, isaresults{i}.seeddata);
  N = N + isaresults{i}.rundata.N;
end
isaresults{1}.rundata.N=N;
result = struct('Rows',rows,'Columns',cols,'seeddata',seeddata,'rundata',isaresults{1}.rundata); 


% Another filtering 
result = isaUnique(normed_data, result);

% Final Result (biClustResult form)
RowxNumber = bsxfun(@ne, result.Rows, 0);
NumberxCol = bsxfun(@ne, result.Columns, 0)'; % See transposed
Number = size(result.Rows, 2);
clust = struct('rows',[], 'cols',[]);
for i=1:Number
 r = find(RowxNumber(:,i)>0);
 c =find(NumberxCol(i,:)>0);
 clust(i) = struct('rows', r, 'cols', c);
end

biClustResult = struct('RowxNum',RowxNumber, 'NumxCol',NumberxCol, ...
                       'ClusterNo',Number, 'Clust', clust, ...
                       'Rundata', result.rundata, 'Seeddata', result.seeddata);
end
