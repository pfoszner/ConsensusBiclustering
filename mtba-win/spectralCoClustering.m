function biClustResult = spectralCoClustering(input,numClust,display)
%Function to compute clusters of the given matrix using Bipartite Spectral Graph Partitioning
% Idea is to model the data matrix as bipartite graph with two sets of 
% vertices as row and column variables of the matrix. Hence the simultaneous 
% clustering problem reduces to bipartite graph partitioning problem.
% The algorithm follows the paper by:
% Dhillon, Inderjit S. 
% "Co-clustering documents and words using bipartite spectral graph partitioning." 
% In Proceedings of the seventh ACM SIGKDD international conference on 
% Knowledge discovery and data mining, pp. 269-274. ACM, 2001.
% 
% Inputs:
%        input           :    Input matrix whose clusters are to be found. 
%                             There is no restriction on the format of
%                             the matrix but preferably it should be of binary form. 
%        numClust        :    Number of biclusters to be found.
%        display         :    It is an optional parameter to display the 
%                             heatmap showing all biclusters.
%                             Value : 1 (if heatmap is required)  
%                                   : 0 (if heatmap is not required) (Default value)
%                               
% Outputs:
%         biClustResult :   A structure consisting of
%           RowxNum   - Structure consisting of Logical Matrix for positive and negative data.
%                       Value is 1 in [i,j] if Row i is in Bicluster j.
%           NumxCol   - Structure consisting of Logical Matrix for positive and negative data. 
%                       Value is 1 in [i,j] if Col j is in Bicluster i.
%           ClusterNo - It is of the form; [Number of Bi-clusters in Positive data, Number of Bi-clusters in Negative data]
%           Clust     - Another structure array containing all clusters with
%                       their respective row and column indices.
% 
% See Also: SpectralHeatmap
%
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Input argument check
if nargin < 1
    error('No Input argument specified.');
end
if nargin < 2
    error('Number of Biclusters not specified.');
end
if nargin < 3
    display = 0;
end
%% Executing the algorithm
matrix=binarize(input); % Binarizing the input matrix
D1=diag(sum(matrix,2)); % Sum each row
D2=diag(sum(matrix,1)); % Sum each column

%Approximate zero values by epsilon
D1(find(D1==0))=eps;   
D2(find(D2==0))=eps;

%Step1: Form the A_n (normalised matrix) 
D1_root=abs((D1)^(-0.5));
D2_root=abs((D2)^(-0.5));
An=D1_root*matrix*D2_root; 

%Step2: Get log2(n) singular vectors and form the Z matrix
[U,S,V]=svd(An,'econ'); 
num_of_pcs=ceil(log2(numClust));
if numClust==1
    num_of_pcs=1;
end
z=[D1_root*U(:,2:num_of_pcs+1);D2_root*V(:,2:num_of_pcs+1)]; 

%Step3: Run k-means
clust_idx=kmeans(z,numClust); 
row_cluster=clust_idx(1:size(matrix,1));
column_cluster=clust_idx(size(matrix,1)+1:end);

rowxNum = zeros(size(matrix,1),numClust);
numxCol = zeros(numClust,size(matrix,2));
for i=1:size(rowxNum,1)
    rowxNum(i,row_cluster(i))=1;
end
for i=1:size(numxCol,2)
    numxCol(column_cluster(i),i)=1;
end

cluster = struct('rows', [], 'cols',[]);
for i=1:numClust
    rows = find(rowxNum(:,i)>0);
    cols = find(numxCol(i,:)>0);
    cluster(i) = struct('rows', rows, 'cols', cols);
end

biClustResult = struct('RowxNum', rowxNum, 'NumxCol', numxCol, 'ClusterNo', numClust, ...
                'Clust', cluster);
%If display is 1 then show the heat map
if display
    SpectralHeatmap(input,row_cluster,column_cluster);
end
end