function biClustResult = LAS(matrix, numBC, threads, iterationsPerBC, scoreThreshold)
%Function to compute bi-clusters using LAS(Large Average Submatrices) method
%  LAS is a statistically motivated biclustering algorithm that tries to
%  find large average submatrices within the data matrix. A score function
%  based on Bonferroni significance correction is also defined which trades
%  off between the submatrix size and average value.
%  The algorithm follows the paper by:
%  Shabalin, Andrey A., Victor J. Weigman, Charles M. Perou, and Andrew B.
%  Nobel. "Finding large average submatrices in high dimensional data." 
%  The Annals of Applied Statistics (2009): 985-1012.
%  
% Inputs: 
%   matrix          = contains data whose bi-clusters are to be found
%   numBC           = Number of Bi-clusters required (aprioi information)
%   iterationsPerBC = Number of iterations per Bi-cluster
%   scoreThreshold  = Threshold value of score. If the value of score
%                     s below the threshold value then the iterations are 
%                     stopped. 
%   threads         = option for parallel computing in MATLAB
%                     Set the value greater than 1.
% Outputs: 
%   biClustResult  : A structure consisting of
%           RowxNum   - Structure consisting of Logical Matrix for positive 
%                       and negative data.
%                       Value is 1 in [i,j] if Row i is in Bicluster j.
%           NumxCol   - Structure consisting of Logical Matrix for positive
%                       and negative data. 
%                       Value is 1 in [i,j] if Col j is in Bicluster i.
%           ClusterNo - It is of the form; [Number of Bi-clusters in 
%                       Positive data, Number of Bi-clusters in Negative data]
%           Clust     - Another structure array containing all clusters with
%                       their respective row and column indices.
%
% See Also: LAS_score
%
% Author : Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Input argument check
if nargin < 1
    error('No Input specified.');
end
if nargin < 2
    numBC = 100;
    iterationsPerBC = 1000;
    scoreThreshold = 1;
    threads = 1;
end 
if nargin < 3
    iterationsPerBC = 1000;
    scoreThreshold = 1;
    threads = 1;
end
if nargin < 4
    iterationsPerBC = 1000;
    scoreThreshold = 1;
end
%% Main routine for LAS
if(threads>1)
	if(matlabpool('size')==0)
		matlabpool 1
	end;
end;
    
BCnumPositive = numBC;  BCnumNegative = numBC;
[rowClust, columnClust, numBiclusters] = LASmainfile(matrix, BCnumPositive, BCnumNegative, iterationsPerBC, scoreThreshold);
RowxNum = rowClust(1:numBiclusters,:)';
ColxNum = columnClust(1:numBiclusters,:);
clusterNo = numBiclusters;

for j=1:size(RowxNum,2)
    rowsPos = find(RowxNum(:,j)>0);
    colsPos =find(ColxNum(j,:)>0);
    cluster(j) = struct('rows', rowsPos, 'cols', colsPos); 
end

biClustResult = struct('RowxNum',RowxNum, 'NumxCol',ColxNum, 'ClusterNo',clusterNo, 'Clust',cluster);

end
