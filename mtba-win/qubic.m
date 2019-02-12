function biClustResult = qubic(datamatrix, quantile, rank, consistancy, output_bicluster, filter_proportion)
% QUalitative BIClustering Algorithm
%
% This algorithm uses statistical model to find all statistically
% significant biclusters including Biclusters with scalling patterns.
% It follows the paper by :
% Guojun Li, Qin Ma, Haibao Tang, Andrew H. Paterson and Ying Xu 
% "QUBIC: a qualitative biclustering algorithm for analyses of gene
% expression data" Nucleic Acids Research, vol. 37, no. 15, 2009 
%
% Usage
% >> biClustResult =  qubic(datamatrix, quantile, rank, consistancy, output_bicluster, filter_proportion)
%
% Inputs:
%   datamatrix                   - Input matrix
%   quantile(optional)           - Percentage of regulating condition for each
%                                  gene
%   rank(optional)               - Range of possible ranks
%   consistency(optional)        - consistency level
%   output_bicluster(optional)   - desired number of output biclusters
%   filter_proportion(optional)  - central parameter for overlap among
%                                  to-be-defined bicluster
%
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] if Row i is 
%                     in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] if Col j is 
%                     in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices.
%
% See Also :  qubicAlgorithm.m
%  
% Author: Shruti jain, 2014
%        
% Contact: sjain@iitk.ac.in, sweetushruti963@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India
%
%% Input argument check
if nargin<1
    error('No input argument specified');
end
if nargin<2
    quantile = 0.06;
    rank = 1;
    consistancy = 0.95;
    output_bicluster = 100;
    filter_proportion = 1;
end
if nargin<3
    rank = 1;
    consistancy = 0.95;
    output_bicluster = 100;
    filter_proportion = 1;
end
if nargin<4
    consistancy = 0.95;
    output_bicluster = 100;
    filter_proportion = 1;
end
if nargin<5
    output_bicluster = 100;
    filter_proportion = 1;
end
if nargin<6
    filter_proportion = 1;
end
  
%% Executing the main program
model = qubicAlgorithm(datamatrix, quantile, rank, consistancy, output_bicluster, filter_proportion);
RowxNum = zeros(size(datamatrix,1),size(model,2));
NumxCol = zeros(size(model,2),size(datamatrix,2));
ClusterNo = size(model,2);

for i=1:size(model,2)
    RowxNum(logical(model(i).GeneIndex),i) = 1;
    NumxCol(i,logical(model(i).TissueIndex)) = 1;
    Clust(i) = struct('rows',find(model(i).GeneIndex),'cols',find(model(i).TissueIndex));
end
biClustResult = struct('RowxNum',RowxNum,'NumxCol',NumxCol,'ClusterNo',ClusterNo,'Clust',Clust);
end
            

    
    