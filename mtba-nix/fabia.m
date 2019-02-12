function biClustResult = fabia(matrix,bicluster_no,sparseness_factor,iteration_no,spl,spz,non_negative,random,center,norm,scale,lap,nL,lL,bL)

% Function to compute biclusters using fabia(Factor Analysis for Bicluster Acquistion)
%
% This algorithm follows the paper by :
% Hochreiter, Sepp, et al. "FABIA: factor analysis for bicluster acquisition." 
% Bioinformatics 26.12 (2010): 1520-1527.
%
%  Usage >> biClustResult = fabia(matrix,bicluster_no,sparseness_factor,iteration_no)
% 
% Inputs :
%             matrix  -  Input datamatrix
%       bicluster_no  -  Number of Biclusters to be found;default = 5
%  sparseness_factor  -  sparseness loadings;default = 0.1
%       iteration_no  -  Number of Iterations;default = 500
%                spl  -  sparseness prior loadings;default = 0;
%                spz  -  sparseness factor;default = 0.5
%       non_negative  -  non_negative loading and factors :0 for no, 1 for yes:default = 0
%             random  -  initialization of loadings : <=0 for svd, >0 for random initialization;default = 1
%             center  -  data centering : 1 for mean, 2 for median, 3 for mode and 0 for no centering;default = 2
%               norm  -  data normalization : 1 for (0.75-0.25 quantile), >1 for (var=1), 0 for no normalization;default = 1
%              scale  -  loading vectors are scaled in each iteration to the given variance. 0.0 indicates non scaling; default = 0        
%                lap  -  minimal value of the variational parameter, default = 1
%                 nL  -  maximal number of biclusters at which a row element can participate; default = 0(no limit)
%                 lL  -  maximal number of row elements per bicluster; default = 0 (no limit)
%                 bL  -  cycle at which the nL or lL maximum starts; default = 0 (start at the beginning)
%               
% Output : 
%         biClustResult :   A structure consisting of
%           RowxNum   - Structure consisting of Logical Matrix for positive and negative data.
%                       Value is 1 in [i,j] if Row i is in Bicluster j.
%           NumxCol   - Structure consisting of Logical Matrix for positive and negative data. 
%                       Value is 1 in [i,j] if Col j is in Bicluster i.
%        NumxColopp   - Structure consisting of Logical Matrix for positive and negative data. 
%                       Value is 1 in [i,j] if Col j is in Opposite Bicluster i.
%           ClusterNo - It is of the form; [Number of Bi-clusters in Positive data, Number of Bi-clusters in Negative data]
%           Clust     - Another structure array containing all clusters with
%                       their respective row and column indices.
% See Also : preprocessing.m, fabic.m, infocontent.m, extractbic.m
% 
% Author : Shruti Jain, 2014
%
% Contact : sjain@iitk.ac.in, sweetushruti963@gmail.com
%           Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Input argument check
if nargin<1
    error('No input argument specified');
end
if nargin<15
    bL = 0;
end
if nargin<14
    lL = 0;
end
if nargin<13
    nL = 0;
end
if nargin<12
    lap = 1;
end
if nargin<11
    scale = 0;
end
if nargin<10
    norm = 1;
end
if nargin<9
    center = 2;
end
if nargin<8
    random = 1;
end
if nargin<7
    non_negative = 0;
end
if nargin<6
    spz = 0.5;
end
if nargin<5
    spl = 0;
end
if nargin<4
    iteration_no = 500;
end
if nargin<3
    sparseness_factor = 0.1;
end
if nargin<2
    bicluster_no = 5;
end

%% Executing the main algorithm
model = preprocessing(matrix,bicluster_no,sparseness_factor,iteration_no,spl,spz,non_negative,random,center,norm,scale,lap,nL,lL,bL);
rowxnum = zeros(size(matrix,1),size(model,2));
numxcol = zeros(size(model,2),size(matrix,2));
numxcolopp = zeros(size(model,2),size(matrix,2));
clusterno = size(model,2);

for i=1:size(model,2)
    rowxnum(model(i).GeneIndex,i) = 1;
    numxcol(i,model(i).TissueIndexp) = 1;
    numxcolopp(i,model(i).TissueIndexn) = 1;
    clust(i) = struct('rows',model(i).GeneIndex,'cols',model(i).TissueIndexp,'colopp',model(i).TissueIndexn);
end
biClustResult = struct('RowxNum',rowxnum,'NumxCol',numxcol,'NumxColopp',numxcolopp,'ClusterNo',clusterno,'Clust',clust);

end
