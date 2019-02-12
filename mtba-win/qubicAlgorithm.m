function finModel = qubicAlgorithm(datamatrix, quantile, rank, consistancy, output_bicluster, filter_proportion)
%% Function to find biclusters using QUBIC(quantitative biclustering) Algorithm
%
% Usage >> finModel = qubicAlgorithm(datamatrix, quantile, rank, consistancy, output_bicluster, filter_proportion)
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
%         finModel               - Final bicluster model containing the gene
%         
% Author: Shruti jain, 2014
%        
% Contact: sjain@iitk.ac.in, sweetushruti963@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Executing the qubic algorithm

matrix = datamatrix;
nRow = size(matrix,1);
nCol = size(matrix,2);

%% quantile discretization of the data
[mat,index] = sort(matrix,2);
mid = round(nCol/2);
s = round((nCol*quantile)+1);
lower = zeros(1,nRow);
upper = zeros(1,nRow);
for i=1:nRow
    midElement= mat(i,mid);
    vsf = mat(i,s);
    vse = mat(i,nCol-s+1);
    diff = min(midElement-vsf,vse-midElement);count = 0;
    for j=1:nCol
        if matrix(i,j) >= (midElement+diff)
            upper(i) = upper(i)+1;
        else
            if matrix(i,j) <= (midElement-diff)
                lower(i) = lower(i)+1;
            else
                matrix(i,j) = 0;count = count+1;
            end
        end
    end
end

for i=1:nRow
    rprev = 1;
    for j=1:rank
        r = round(lower(i)*j/rank);
        for k=rprev:r
            matrix(i,index(i,k)) = -rank+j-1;
        end
        rprev = r+1;
    end
end
for i = 1:nRow
    rprev = nCol;
    for j = 1:rank
        r = nCol-round(upper(i)*j/rank)+1;
        for k = r:rprev
            matrix(i,index(i,k)) = rank-j+1;
        end
        rprev = r-1;
    end
end

%% finding the weighted graph to seed data
weight = zeros(3,nRow*(nRow-1)/2);
l=1;
for i =1:nRow-1
    for j = i+1:nRow
        temp = compare_row(matrix(j,:),matrix(i,:),(1:nCol));
        weight(1,l) = sum(temp.res);
        weight(2,l) = i;
        weight(3,l) = j;
        l = l+1;
    end
end
weight = sortrows(weight', 1)';
weight = fliplr(weight);
%weight is the sorted matrix in descending order with second row and third
%row represnting the rows which are compared and first row as weight
%of these rows 

%% finding the qubic biclusters

biclust = 0;       
o = output_bicluster;
c = consistancy;
mincol = max(2,floor(0.05*nCol));
model = struct('GeneIndex',[],'TissueIndex',[]);
model.GeneIndex = zeros(1,nRow);
model.TissueIndex = zeros(1,nCol);
seed = nRow*(nRow-1)/2;

for i = 1:seed
    if weight(1,i)<mincol
        break;
    end
    index1 = weight(2,i);
    index2 = weight(3,i);
    arr1 = zeros(1,biclust);
    arr2 = zeros(1,biclust);
    
    for j = 1:biclust
        arr = find(model(j).GeneIndex);
        if find(arr==index1)
            arr1(j) = j;
        end
        if find(arr==index2)
            arr2(j) = j;
        end
    end
    
    if sum(arr1)~=0 && sum(arr2)~=0
        temp = compare_row(arr1,arr2,(1:biclust));
        temp = sum(temp.res);
        if temp~=0
            continue;
        else
            b1 = find(arr1,1);
            b2 = find(arr2,1);
            size_intersect = size(intersect(find(model(b1).GeneIndex),find(model(b2).GeneIndex)),2);
            if size_intersect > 0 
                continue;
            end
            if weight(1,i) < max(nnz(model(b1).TissueIndex),nnz(model(b2).TissueIndex))
                continue;
            end
        end
    end
    
 %% step-2 : Expansion while maintaing total columnwise consistancy
   biclust = biclust+1;
   model(biclust).GeneIndex = zeros(1,nRow);
   model(biclust).TissueIndex = zeros(1,nCol);
   model(biclust).GeneIndex(index1) = index1;
   model(biclust).GeneIndex(index2) = index2;
   row_no = 2;
   col_no = 0;
   for j=1:nCol
       if (matrix(index1,j)==matrix(index2,j))&&(matrix(index1,j)~=0)
           model(biclust).TissueIndex(j) = j;
           col_no = col_no+1;
       end
   end
   cand_thr = max(2,floor(mincol*c));
   gene1 = index1;
   consist = zeros(1,nRow);
   cand = zeros(1,nRow);
   cand(index1) = 1;
   cand(index2) = 1;
   for j= 1:nRow
       if j~=index1 && j~=index2
           arr_temp = compare_row(matrix(gene1,:),matrix(j,:),model(biclust).TissueIndex);
           consist(i) = sum(arr_temp.res);
       end
   end
   con = sort(consist,'descend');
   var = min(100,nRow);
   top = con(var);
   for j = 1:nRow
       if consist(j)<top
           cand(j) = 1;
       end
   end
   for j = 3:nRow
       max_cnt = -1;
       max_i = -1;
       for k = 1:nRow
           if cand(k)==0
               cnt = compare_row(matrix(index1,:),matrix(k,:),model(biclust).TissueIndex);
               cnt = sum(cnt.res);
               if cnt < cand_thr
                   cand(k) = 1;
               end
               if cnt > max_cnt
                   max_cnt = cnt;
                   max_i = k;
               end
           end
       end
       if max_cnt < mincol || max_i == -1 
           break;
       else
           if min(j,max_cnt) >= min(row_no,col_no)
               row_no = j;
               col_no = max_cnt;
               model(biclust).GeneIndex(max_i) = max_i;
               col_eq = compare_row(matrix(gene1,:),matrix(max_i,:),model(biclust).TissueIndex);
               col_noteq  = ~col_eq.res;
               model(biclust).TissueIndex(logical(col_noteq)) = 0;
               cand(max_i) = 1;
           end
       end
   end
   
   %% step-3 : Expansion allowing less than total consistancy
   c = consistancy;
   arr_i = find(model(biclust).GeneIndex);
   size_i = size(arr_i,2);
   t = zeros(1,nCol);
   t(logical(model(biclust).TissueIndex)) = matrix(gene1,logical(model(biclust).TissueIndex));
   % adding the colums
   profile = zeros(nCol,2*rank);
   for j = 1:size_i-1
       for k = 1:nCol
           val = matrix(arr_i(j),k);
           if val<0
               val = rank-val;
           end
           if val ==0
               continue;
           end
           profile(k,val) = profile(k,val)+1;
       end
   end
   threshold = ceil(nnz(model(biclust).GeneIndex)*c);
   for j = 1:nCol
       for k = 1:2*rank
           n = profile(j,k);
           if k > rank
               k1 = -k+rank;
           else
               k1 = k;
           end
           if matrix(arr_i(size_i),j) == k1
               n = n+1;
           end
           if n >= threshold
               model(biclust).TissueIndex(j) = j;
               t(j) = k1;
               break;
           end
       end
   end 
  
   % adding the rows
   size_j = nnz(model(biclust).TissueIndex);
   arr_j = model(biclust).TissueIndex;
   arr_i = model(biclust).GeneIndex;
   threshold = floor(size_j*c);
   
   for j = 1:nRow
       if arr_i(j) == j
           continue;
       end
       con = compare_row(t,matrix(j,:),arr_j);
       con = sum(con.res);
       if threshold <= con
           model(biclust).GeneIndex(j) = j;
       end
   end
   
  %% Step-4 : Expansion by adding oppositely regulated genes
  
  % adding negatively regulated genes
  for j = 1:nRow
       if model(biclust).GeneIndex(j) == j
           continue;
       end
       con = compare_row(matrix(index1,:),-matrix(j,:),arr_j);
       con = sum(con.res);
       if threshold <= con
           model(biclust).GeneIndex(j) = j;
       end
  end
  threshold = ceil(nnz(model(biclust).GeneIndex)*c);
  arr_i = find(model(biclust).GeneIndex);
  size_i = size(arr_i,2);
  profile = zeros(nCol,2*rank);
   for j = 1:size_i
       for k = 1:nCol
           val = matrix(arr_i(j),k);
           if val<0
               val = rank-val;
           end
           if val ==0
               continue;
           end
           profile(k,val) = profile(k,val)+1;
       end
   end
   
   for j = 1:nCol
       for k = 1:2*rank
           if profile(j,k) >= threshold
               model(biclust).TissueIndex(j) = j;
           end
       end
   end
  
  if biclust == 2*o
      break;
  end
end

%% Filtering the overlapping biclusters
%num = min(o,biclust);
i = 1;
j = 1;
arr_row = (1:nRow);
arr_col = (1:nCol);
while i <= biclust && j <= o
    cur_row = size(find(model(i).GeneIndex),2);
    cur_col = size(find(model(i).TissueIndex),2);
    flag = 1;
    k = 1;
    while k < j
       temp = compare_row(model(i).GeneIndex,model(k).GeneIndex,arr_row);
       inter_row = sum(temp.res);
       temp = compare_row(model(i).TissueIndex,model(k).TissueIndex,arr_col);
       inter_col = sum(temp.res);
       if inter_row*inter_col > filter_proportion*cur_row*cur_col
           flag = 0;
           break;
       end
       k = k+1;
    end
    if flag == 1
        finModel(j) = struct('GeneIndex',model(i).GeneIndex,'TissueIndex',model(i).TissueIndex);
        j = j+1;        
    end 
    i = i+1;
end   

end


%% Function to evaluate quantile of a row vector  %
function result = quant(row,nCol,f)

index = f*(nCol-1);
lhs = floor(index);
delta = index-lhs;
if delta ==0
   result = (1-delta)*row(lhs+1);
else
    result = (1-delta)*row(lhs+1) + delta*row(lhs+2);
end
end
  