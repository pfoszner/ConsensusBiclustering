function biClustResult = ITL(matrix,numClust,display)
% Function to find biclusters using Information Theoretic learning
% The algorithm was developed for simultaneously clustering the rows and
% columns of contingency table. It views contingency tables as joint
% probability distribution of the two discrete random variables. Hence the
% problem of simultaneously clustering reduces to maximizing the mutual
% information between the clustered random variables or minimizing the loss
% in mutual information as a result of clustering.
% The algorithm follows the paper by:
% Dhillon, Inderjit S., Subramanyam Mallela, and Dharmendra S. Modha.
% "Information-theoretic co-clustering." In Proceedings of the ninth ACM SIGKDD 
% international conference on Knowledge discovery and data mining, pp. 89-98. ACM, 2003.
% 
% Inputs:
%        matrix          =   Input matrix which represents the joint probability
%                            distribution of the two variables
%        numClust        =   Desired number of clusters
%        display         =   Option for plotting the final biclusters on a
%                            heatmap (1 for when required and 0 when not required).
%           
% Output:
%        biClustResult :   A structure consisting of
%           RowxNum   - Structure consisting of Logical Matrix for positive and negative data.
%                       Value is 1 in [i,j] if Row i is in Bicluster j.
%           NumxCol   - Structure consisting of Logical Matrix for positive and negative data. 
%                       Value is 1 in [i,j] if Col j is in Bicluster i.
%           ClusterNo - It is of the form; [Number of Bi-clusters in Positive data, Number of Bi-clusters in Negative data]
%           Clust     - Another structure array containing all clusters with
%                       their respective row and column indices.
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
    numClust = 100;
    display = 0;
end
if nargin < 3
    display = 0;
end
%% Calling the main program
RowClusterNo = numClust;
ColClusterNo = numClust;
[row_cluster,column_cluster] = InformationTheoreticLearning(matrix,RowClusterNo,ColClusterNo);
% making the output structure
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

if display==1
    SpectralHeatmap(matrix, row_cluster,column_cluster);
end
end
