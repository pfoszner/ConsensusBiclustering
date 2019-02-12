function biClustResult = opsm(matrix,partialModelSize)
% Function to compute biclusters using OPSM (Order Preserving Sub-matrices) 
% Ben-Dor et al. defined a bicluster as an order preserving submatrix (OPSM).
% According to their definition, a bicluster is a group of rows whose
% values induce same linear ordering across a subset of the columns. More
% specifically the algorithm aims at finding large OPSMs. A submatrix is
% order preserving if there is a permutation of columns for which the
% values in every row is strictly increasing. For eg:
% Let U be a submatrix = { 3  1  6  7      
%                          4  1  5  8       
%                          7  2  9  11      
%                          5  3  6  7  }   
% Let C_1, C_2, C_3 and C_4 be the four columns, then we can notice that if
% we arrange the columns as C_2, C_1, C_3 and C_4, then we obtain increasing 
% sequences along the rows. Therefore U is an order preserving submatrix.
%
% The algorithm follows the paper by : 
% Ben-Dor, Amir, Benny Chor, Richard Karp, and Zohar Yakhini. "Discovering
% local structure in gene expression data: the order-preserving submatrix problem." 
% In Proceedings of the sixth annual international conference on Computational biology, pp. 49-57. ACM, 2002.
% 
% Inputs :
%         matrix             :   Input matrix [nxm] of co-occurences of n instances (rows) and m
%                                features (columns).
%         partialModelSize   :   Input argument for maximum holding capacity of partial models 
%                                i.e. maximum number of elements allowed in a partial model.
%   
% Output : 
%         biClustResult :   A structure consisting of
%           RowxNum   - Structure consisting of Logical Matrix for positive and negative data.
%                       Value is 1 in [i,j] if Row i is in Bicluster j.
%           NumxCol   - Structure consisting of Logical Matrix for positive and negative data. 
%                       Value is 1 in [i,j] if Col j is in Bicluster i.
%           ClusterNo - It is of the form; [Number of Bi-clusters in Positive data, Number of Bi-clusters in Negative data]
%           Clust     - Another structure array containing all clusters with
%                       their respective row and column indices.
% 
% See Also: opsmAlgorithm
%
% Author : Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Input argument check
if nargin<1
    error('No Input argument defined.');
end
if nargin<2
    partialModelSize = 10;
end
%% Executing the main program
model = opsmAlgorithm(matrix,partialModelSize);
RowxNum = zeros(size(matrix,1),size(model,2));
NumxCol = zeros(size(model,2),size(matrix,2));
ClusterNo = size(model,2);

for i=1:size(model,2)
    RowxNum(model(i).GeneIndex,i) = 1;
    NumxCol(i,model(i).TissueIndex) = 1;
    Clust(i) = struct('rows',model(i).GeneIndex,'cols',model(i).TissueIndex);
end
biClustResult = struct('RowxNum',RowxNum,'NumxCol',NumxCol,'ClusterNo',ClusterNo,'Clust',Clust);
end
