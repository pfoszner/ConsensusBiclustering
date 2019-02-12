function [row_clust,column_clust]=InformationTheoreticLearning(ProbDist,k,l)
%  Implementation of the ITL algorithm
% Step1: Initialization
% Initialize the partition functions for rows and columns
% Step2: Reclustering
% If the difference in the loss of mutual information is not below the
% threshold,then again compute the row and column clusters.
% Inputs:
%        ProbDist      = Matrix consisting of joint probability distribution 
%                        of two discrete random variables.
%           k          = Number of required clusters of the row variable.
%           l          = Number of required clusters of the column variable.
% Output:
%        row_clust     = Vector containing the cluster number of each row.
%        column_clust  = Vector containing the cluster number of each column.
%
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Initialization
 row_clust=[];
 for i=1:k-1
    cl=i*ones(floor(size(ProbDist,1)/k),1);
    row_clust=[row_clust ; cl];
 end
 row_clust=[row_clust ; k*ones(size(ProbDist,1) - (k-1) * floor(size(ProbDist,1)/k),1)];
column_clust=[];
 for i=1:l-1
     cl=i*ones(floor(size(ProbDist,2)/l),1);
     column_clust=[column_clust ; cl];
 end
 column_clust=[column_clust ; l*ones(size(ProbDist,2) - (l-1) * floor(size(ProbDist,2)/l),1)];
%% Reclustering
 loss = 1;
 while loss > 0.01
     init_kl_div = 0;
     fin_kl_div = 0;
     [row_clust,column_clust,q_init,q] = Reclustering(ProbDist,k,l,row_clust,column_clust);      
     for i=1:size(ProbDist,1);
         for j=1:size(ProbDist,2);
             if ProbDist(i,j)~=0
                 init_kl_div = init_kl_div + (ProbDist(i,j) * log2(ProbDist(i,j)/q_init(i,j)));
                 fin_kl_div = fin_kl_div + (ProbDist(i,j) * log2(ProbDist(i,j)/q(i,j)));
             end
         end
     end
     loss = init_kl_div - fin_kl_div;
 end
end
