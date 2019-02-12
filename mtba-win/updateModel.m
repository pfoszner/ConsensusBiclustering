function [model,isInG,isInT] = updateModel(model,score,noCol,lowestRankCol,highestRankCol,capacity,noSmallElements,noLargeElements,extColumn,isInG,isInT,index)
% Function to update the partial model
% The partial model is stored in the form of binary tree. A new element can
% be added to the model in two ways.
% First :
% If the heap is not full, new element is added using putAsUsual function.
% Second : 
% If the heap is full, new element is added only if it has score more than
% the element with lowest score.
% Inputs :
%        model            :    Present Partial Model  
%        score            :    Significance score for new element
%        noCol            :    Number of columns considered
%        lowestRankCol    :    Column with lowest rank
%        highestRankCol   :    Column with highest rank
%        capacity         :    Maximum holding capacity for partial models
% If partial model is of order (a,b), then 
%        noSmallElements  :    a  
%        noLargeElements  :    b
%        extColumn        :    Column being considered for extending the partial model
%        isInG            :    (i,j) entry is 1, if i-th model contains j-th row 
%        isInT            :    (i,j) entry is 1, if i-th model contains j-th column
%        index            :    Partial model number
%        
% Outputs :
%         model           :    Output Partial Model
%         isInG           :    (i,j) entry is 1, if i-th model contains j-th row 
%         isInT           :    (i,j) entry is 1, if i-th model contains j-th column
%        
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, 
%          Indian Institute of Technology Kanpur, India
%% Main routine for updating partial model
% Elements with zero score are not added to the model
if score==0
    model = model;
else
    % If the heap is not full, the new element is added using 'putAsUsual' function
    if size(model,2) < capacity
        [model,isInG,isInT] = putAsUsual(model,score,noCol,lowestRankCol,highestRankCol,noSmallElements,noLargeElements,extColumn,isInG,isInT,index);
        % If the heap is full and the new element has a score value more than that of the 
        % element with lowest score, then the element is added using 'putByReplacingRoot' function
    else if score > model(2).score
        [model,isInG,isInT] = putByReplacingRoot(model,score,noCol,lowestRankCol,highestRankCol,noSmallElements,noLargeElements,extColumn,isInG,isInT,index);
        end
    end
end
end

function [mod,isInG,isInT] = putAsUsual(mod,p,s,ta,tb,x,y,t,isInG,isInT,index)
%% Function to put new element to the heap when the heap is not full
% This function is used to add new element to the model if the heap is not full.
%  
% Inputs :
%         mod       :    Present partial model
%         p         :    Significance score for new element
%         s         :    Number of columns considered
%         ta        :    Column with lowest rank
%         tb        :    Column with highest rank
% If partial model is of order (a,b), then 
%         x         :    a  
%         y         :    b
%         t         :    Column being considered for extending the partial model
%         isInG     :    (i,j) entry is 1, if i-th model contains j-th row
%         isInT     :    (i,j) entry is 1, if i-th model contains j-th column
%         index     :    Partial model number
% Outputs :
%         mod       :    Ouput partial model
%         isInG     :    (i,j) entry is 1, if i-th model contains j-th row 
%         isInT     :    (i,j) entry is 1, if i-th model contains j-th column
%% Main routine of the function
% Checking whether present model is empty or not.
if isempty(mod)
        element = struct('score',p,'noCol',s,'lowestRankCol',ta,'highestRankCol',tb,'noSmallElements',x,'noLargeElements',y,'extColumn',t,'modelNo',index);
        mod = struct('score',[],'noCol',[],'lowestRankCol',[],'highestRankCol',[],'noSmallElements',[],'noLargeElements',[],'extColumn',[],'modelNo',[]);
        mod = [mod,element];
% If the present model is not empty, then add the element to it's correct
% position in the heap.
else
    element = struct('score',p,'noCol',s,'lowestRankCol',ta,'highestRankCol',tb,'noSmallElements',x,'noLargeElements',y,'extColumn',t,'modelNo',index);
    % Initially adding the element to the end of the model
    mod = [mod,element];
    len = size(mod,2);
    i = len;
    parent = floor((i+1)/2); 
    term = mod(i).score < mod(parent).score;
    if isempty(term)
        term = 0;
    end
    % Let the new element dive up to it's correct position
    while ((parent>=2) && term)
        temp = mod(parent);
        mod(parent) = mod(i);
        mod(i) = temp;        
        i = parent;
        parent = floor((i+1)/2);
        term = mod(i).score < mod(parent).score;
        if isempty(term)
            term = 0;
        end
    end
    clear term;
end
end

function [mod,isInG,isInT] = putByReplacingRoot(mod,p,s,ta,tb,x,y,t,isInG,isInT,index)
%% Function to put new element to the heap when the heap is full
% This function removes the element with lowest score from the model and adds the new element.
%
%  Inputs : 
%         mod       :    Present partial model
%         p         :    Significance score for new element
%         s         :    Number of columns considered
%         ta        :    Column with lowest rank
%         tb        :    Column with highest rank
% If partial model is of order (a,b), then 
%         x         :    a  
%         y         :    b
%         t         :    Column being considered for extending the partial model
%         isInG     :    (i,j) entry is 1, if i-th model contains j-th row
%         isInT     :    (i,j) entry is 1, if i-th model contains j-th column
%         index     :    Partial model number
% Outputs :
%         mod       :    Ouput partial model
%         isInG     :    (i,j) entry is 1, if i-th model contains j-th row 
%         isInT     :    (i,j) entry is 1, if i-th model contains j-th column
%% Main routine of the function
newElement = struct('score',p,'noCol',s,'lowestRankCol',ta,'highestRankCol',tb,'noSmallElements',x,'noLargeElements',y,'extColumn',t,'modelNo',index);
% Replacing the root of the model with new element
mod(2) = newElement;
i = 2;
right = (i*2); 
left =(i*2)-1;

if right > size(mod,2)
    ter = 0;
else
    ter = mod(i).score > mod(right).score;
end
% Let new element sink down either left or right as long as both entries exist
while ((right<=size(mod,2)) && mod(i).score>mod(left).score || ter)
    if mod(left).score>mod(right).score
        temp = mod(right);
        mod(right) = mod(i);
        mod(i) = temp;
        i = right;
    else
        temp = mod(left);
        mod(left) = mod(i);
        mod(i) = temp;
        i = left;
    end
    right = (i*2);
    left = (i*2)-1;
    if right > size(mod,2)
        ter = 0;
    else
        ter = mod(i).score > mod(right).score;
    end
end
% Let new element sink down left when only left entries exist
while (left<=size(mod,2) && mod(i).score>mod(left).score)   
    temp = mod(left);
    mod(left) = mod(i);
    mod(i) = temp;
    i = left;
    left = (i*2) - 1;
end
end
