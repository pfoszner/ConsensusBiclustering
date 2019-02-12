function [row_clust,column_clust,q_init,q] = Reclustering(ProbDist,RowClusterNo,ColClusterNo,row_clust,column_clust)
% Function to perform row and column reclustering
% Step1: Compute all the required probabilities
% The probabilities required for row reclustering are p(Y|x) and q(Y|xhat)
%
% Step2: Perform row reclustering
% Assign new cluster values to each row, which minimizes the KL divergence
% of p(Y|x) and q(Y|xhat).
% 
% Step3: Compute the probabilities according to new row cluster values 
% The probabilities required for column reclustering are p(X|y) and q(X|yhat)
% 
% Step4: Perform column reclustering
% Assign new cluster values to each column, which minimizes the KL divergence
% of p(X|y) and q(X|yhat).
%
% Inputs :
%        ProbDist = Input Probability Distribution matrix; p(X,Y)
%        RowClusterNo   =   Number of clusters required for row variable       
%        ColClusterNo   =   Number of clusters required for column variable
%        row_clust      =   Vector containing cluster number for each row
%        column_clust   =   Vector containing cluster number for each column
% Output :
%        row_clust      =   Vector containing cluster number for each row
%        column_clust   =   Vector containing cluster number for each column
%        q_init         =   q(X,Y) before performing reclustering
%           q           =   q(X,Y) obtained after performing row and column reclustering
% 
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Step1: Compute all the required probabilities for row reclustering
k = RowClusterNo;
l = ColClusterNo;
% computing q(X,Y)
q = updateQ(ProbDist,row_clust,column_clust,k,l);
q_init = q;

% computing q(X_hat,Y_hat)
q_hat = zeros(k,l);
for i = 1:k
    for j = 1:l
        q_hat(i,j) = sum(sum(q(find(row_clust == i), find(column_clust == j))));
    end
end

% computing q(Y|X_hat) = q(X,Y)/q(X_hat)
q_yGivXhat=zeros(k,size(ProbDist,2));
for i=1:k
    for j=1:size(q_yGivXhat,2)
        ind = find(row_clust==i);
        q_yGivXhat(i,j) = sum(sum(q(ind,j))) / sum(sum(ProbDist(ind,:)));
    end
end

%% Step2 : Perform Row Re-clustering
% Re-assigning the cluster values to each row, according to following expression
%    row_clust(i) = Value of X_hat that minimizes (D(p(Y|X) || q(Y|X_hat)))
for i=1:size(ProbDist,1)
    for j=1:k
        arg(i,j)=0;
        for temp=1:size(ProbDist,2)
            if ProbDist(i,temp)==0
                arg(i,j)=arg(i,j);
            else
                arg(i,j) = arg(i,j) + (ProbDist(i,temp) / sum(ProbDist(i,:))) * log2(ProbDist(i,temp) / (q_yGivXhat(j,temp) * sum(ProbDist(i,:))));
            end
        end
    end
    row_clust(i)=min(abs(arg(i,:)));
    if length(find(abs(arg(i,:))==row_clust(i))) > 1
        temp = find(abs(arg(i,:))==row_clust(i));
        row_clust(i) = temp(1);
    else
        row_clust(i)=find(abs(arg(i,:))==row_clust(i));
    end
    clear temp;
end

%% Step3 : Compute the probabilities according to new row cluster values
% Updating q(X,Y)
q = updateQ(ProbDist,row_clust,column_clust,k,l);
% Updating q(X_hat,Y_hat)
q_hat = zeros(k,l);
for i = 1:k
    for j = 1:l
        q_hat(i,j) = sum(sum(q(find(row_clust == i), find(column_clust == j))));
    end
end
% Computing q(X|Y_hat) = q(X,Y)/q(Y_hat)
q_xGivYhat = zeros( size(ProbDist,1) , l);
for i=1:size(q_xGivYhat,1)
    for j=1:l
        ind_1 = find(column_clust==j);
        q_xGivYhat(i,j) = sum(sum(q(i,ind_1)))/sum(sum(ProbDist(:,ind_1)));
    end
end
clear arg;

%% Step4 : Column Re-clustering
% Re-assigning the cluster values to each column, according to following expression
%    column_clust(i) = Value of Y_hat that minimizes (D(p(X|Y) || q(X|Y_hat))) 
for i = 1:size(ProbDist,2)
    for j = 1:l
        arg(i,j)=0;
        for temp = 1:size(ProbDist,1)
            if ProbDist(temp,i) == 0
                arg(i,j)=arg(i,j);
            else
                arg(i,j) = arg(i,j) + (ProbDist(temp,i) / sum(ProbDist(:,i))) * log2(ProbDist(temp,i) / (q_xGivYhat(temp,j) * sum(ProbDist(:,i))));
            end
        end
    end
    column_clust(i)=min(abs(arg(i,:)));
    if length(find(abs(arg(i,:))==column_clust(i))) > 1
        temp = find(abs(arg(i,:))==column_clust(i));
        column_clust(i) = temp(1);
    else
        column_clust(i)=find(abs(arg(i,:))==column_clust(i));
    end
    clear temp;
end

% Updating q(X,Y)
q = updateQ(ProbDist,row_clust,column_clust,k,l);
end

function q = updateQ(ProbDist,row_clust,column_clust,RowClusterNo,ColClusterNo)
%% Function to compute q(X,Y) using p(X,Y), row_clust and column_clust
% q(X,Y) is computed using following expression:
%     q(X,Y) = p(X_hat,Y_hat)*p(X|X_hat)*p(Y|Y_hat)
%
% Inputs :
%        ProbDist     =   Input Probability Distribution; p(X,Y)
%        row_clust    =   Vector containing cluster number for each row
%        column_clust =   Vector containing cluster number for each column
%        RowClusterNo =   Number of clusters required for rows
%        ColClusterNo =   Number of clusters required for columns
% Output : 
%        q            =   q(X,Y)
%% Main routine for 'updateQ'
k = RowClusterNo;
l = ColClusterNo;
% Computing p(X_hat,Y_hat)
ProbDist_hat = zeros(k,l);
for i = 1:k
    for j = 1:l
        ProbDist_hat(i,j) = sum(sum(ProbDist(find(row_clust == i), find(column_clust == j)))); 
    end
end
% Computing p(X|X_hat) = p(X)/p(X_hat)
p_xGivXhat = zeros(size(ProbDist,1),k);
for i=1:size(ProbDist,1)
    p_xGivXhat(i,row_clust(i))=sum(ProbDist(i,:))/sum(ProbDist_hat(row_clust(i),:));
end
% Computing p(Y|Y_hat) = p(Y)/P(Y_hat)
p_yGivYhat = zeros(l,size(ProbDist,2));
for i=1:size(ProbDist,2)
    p_yGivYhat(column_clust(i),i) = sum(ProbDist(:,i))/sum(ProbDist_hat(:,column_clust(i)));
end
% Computing q(X,Y)
q = p_xGivXhat * ProbDist_hat * p_yGivYhat;
end
