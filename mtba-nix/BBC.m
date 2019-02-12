function biClustResult = BBC(matrix,numClusters,normChoice,alpha)
% Function to find biclusters using Bayesian Biclustering
% This technique uses Gibbs Sampling and Bayesian Information Criterion to
% identify biclusters. Gibbs Sampling helps in drawing inference from BBC
% model. The algorithm follows the paper:
% Gu, Jiajun, and Jun S. Liu. "Bayesian biclustering of gene expression
% data." BMC genomics 9, no. Suppl 1 (2008): S4.
% Inputs:
%        matrix         :     Input dataset whose biclusters are to be found
%        numClusters    :     Number of biclusters required
%        normChoice     :     Parameter for selecting normalization method
%                             0 : No normalization required (But it is recommended to use some normalization technique)
%                             1 : Row Standardization (default value)
%                             2 : Column Standardization
%                             3 : Interquartile Range Normalization
%                             4 : Smallest Quartile Range Normalization
%        alpha          :     It is an optional parameter and required only
%                             when normChoice is either 3 or 4. Value should 
%                             be between 5 and 100.
% Outputs:
%   biClustResult  : A structure consisting of following fields
%           RowxNum   - Logical Matrix 
%                       Value is 1 in [i,j] if Row i is in Bicluster j.
%           NumxCol   - Logical Matrix 
%                       Value is 1 in [i,j] if Col j is in Bicluster i.
%           ClusterNo - Number of Bi-clusters found.
%           Clust     - Another structure array containing all clusters with
%                       their respective row and column indices.
%
% Author : Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

  %% Input argument check
if nargin<1
    error('All input arguments not defined.');
end
if nargin<2
    numClusters=100;
    normChoice=1;
    alpha=0;
end
if nargin<3
    normChoice=1;
    alpha=0;
end
if nargin<4
    if normChoice==3 || normChoice==4
        alpha=90;
    else
        alpha=0;
    end
end
%% Calling the main program
[rowxNum, colxNum] = mBBC(matrix,numClusters,normChoice,alpha);
for i=1:numClusters
    rows = find(rowxNum(i,:));
    col = find(colxNum(i,:));
    clust(i)= struct('rows',rows,'cols',col);
end
biClustResult = struct('RowxNum',rowxNum','NumxCol',colxNum,'ClusterNo',numClusters,'Clust',clust);
end
