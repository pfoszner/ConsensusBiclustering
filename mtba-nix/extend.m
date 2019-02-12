function [isInT,arrayTissue,parent,isInG] = extend(dataset,extColumn,model,arrayTissue,isInT,isInG)
% Function to realize best l models
% This function realizes the best l models found using significance score.
% The model is saved in 'parent' which consists of appropriate fields.
%
% Inputs: 
%        dataset          :   Input matrix
%        extColumn        :   Column being considered for extending the partial model
%        model            :   Model to be realized
%        arrayTissue      :   Array containing column numbers present in the models 
%        isInT            :   (i,j) entry is 1, if i-th model contains j-th column
%        isInG            :   (i,j) entry is 1, if i-th model contains j-th row 
% Outputs:       
%        isInT            :   (i,j) entry is 1, if i-th model contains j-th column
%        arrayTissue      :   Array containing column numbers present in the models 
%        parent           :   Array of best l realized models 
%        isInG            :   (i,j) entry is 1, if i-th model contains j-th row 
%
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of
%          Technology, Kanpur, India
%% Main routine for the function
t = extColumn;
% The following if-else condition is for checking whether the model is being realized for first time or not 
if isempty(arrayTissue)
    [parent,isInT,arrayTissue,isInG] = realize(dataset,model);
else
    noRows = size(dataset,1);
    a=model.noSmallElements;
    b=model.noLargeElements;

    if a+b==length(arrayTissue)
        error('Complete Models cannot be extended!');
    end
    % Following if-else condition appropriately increases 'a' or 'b' parameters of the partial model
    if a==b
        arrayTissue(a+1) = t; a=a+1;
    else
        arrayTissue(length(arrayTissue)-b) = t; b = b+1;
    end
    [so,D] = sort(dataset,2);
    [s1,ind] = sort(D,2);
    isInT(t) = 1;
    % Checking compatibility of each row with the model. isInG[] is updated 
    isInG = updateCompatibilityInfo(ind,noRows,arrayTissue,a,b,isInG);
    model.noSmallElements = a; model.noLargeElements = b;
    parent = model;
end
end

function [finalmodel,isInT,arrayTissue,isInG] = realize(dataset,element)
%% Function to realize the initial partial model before extending
% This function is used for realizing the model which has not been extended yet.
% Inputs:
%         dataset      :     Input dataset
%         element      :     Input partial model to be realized
% Outputs:
%         finalmodel   :     Final model obtained 
%         isInT        :     (i,j) entry is 1, if i-th model contains j-th column
%         arrayTissue  :     Array containing column numbers present in the models 
%         isInG        :     (i,j) entry is 1, if i-th model contains j-th row 
%% Main routine for the function
% Computing 'arrayTissue' vector
arrayTissue(1) = element.lowestRankCol;
arrayTissue(element.noCol) = element.highestRankCol;
arrayTissue(2:(element.noCol-1)) = -1*ones(1,element.noCol-2);
a=1;b=1;
% Computing 'isInT' vector
isInT = zeros(1,size(dataset,2));
isInT(element.lowestRankCol)=1;isInT(element.highestRankCol)=1;
[so,D] = sort(dataset,2);
[s1,ind] = sort(D,2);
% Computing 'isInG' vector
isInG = zeros(1,size(dataset,1));
for i=1:size(dataset,1)
    isInG(i) = ind(i,element.lowestRankCol)<ind(i,element.highestRankCol);
end
% Updating 'isInG'
isInG = updateCompatibilityInfo(ind,size(dataset,1),arrayTissue,a,b,isInG);
arrayGeneIndex = updateG(size(dataset,1),isInG);
finalmodel = struct('geneIndex',arrayGeneIndex,'tissueIndex',find(isInT),'noSmallElements',a,'noLargeElements',b);
end

function isInG = updateCompatibilityInfo(ind,noRows,arrayTissue,a,b,isInG)
%% Function to check for each row if it is compatible with the model
% Inputs:
%         ind             :   Ranked input matrix i.e (i,j) entry of 'ind' corresponds
%                             to the position of (i,j) entry of input dataset 
%                             after sorting the matrix along rows
%         noRows          :   Number of rows in the input matrix
%         arrayTissue     :   Array containing column numbers present in the models 
% If partial model is of order (x,y), then 
%         a               :   Number of smaller elements i.e x 
%         b               :   Number of larger elements i.e y
%         isInG           :   (i,j) entry is 1, if i-th model contains j-th row
% Outputs:       
%         isInG           :   (i,j) entry is 1, if i-th model contains j-th row
%% Main routine of the function
lengthArrayGene = 0;
for i=1:noRows
    if isInG(i)
        isInG(i) = (a==1 || ind(i,arrayTissue(a-1))<ind(i,arrayTissue(a))) ...
            && (ind(i,arrayTissue(a))<ind(i,arrayTissue(length(arrayTissue)-b+1)))...
            && (b==1 || ind(i,arrayTissue(length(arrayTissue)-b+1))<ind(i,arrayTissue(length(arrayTissue)-b+2)));
        if isInG(i)
            lengthArrayGene = lengthArrayGene + 1;
        end
    end
end
end

function arrayGeneIndex = updateG(noRows,isInG)
%% Function to compute array of Genes present in the model according to updated 'isInG'
% This function basically stores the gene indices present in the model.
% Inputs:
%         noRows           :    Number of rows present in the input dataset
%         isInG            :    (i,j) entry is 1, if i-th model contains j-th row
% Outputs:
%         arrayGeneIndex   :    Array containing gene indices present in the models 
%% Main routine of the function
cursor = 1;
arrGene = [];
for i=1:noRows
    if isInG(i)
        arrGene(cursor) = i;
        cursor = cursor + 1;
    end
end
arrayGeneIndex = arrGene;
end
