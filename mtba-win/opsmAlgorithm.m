function finModel = opsmAlgorithm(matrix,partialModelSize)
% Function to find biclusters using OPSM (Order Preserving Sub-Matrices) Algorithm
% Inputs:
%        matrix           :       Input matrix [n*m] of co-occurences of n instances (rows) and m
%                                 features (columns).
%        partialModelSize :       Input argument for maximum capacity of partial models 
%                                 i.e. maximum number of elements allowed in a partial model.
%   
% Outputs:
%         finModel        :       Final bicluster model containing the gene
%         
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India
%% Main routine for 'opsmAlgorithm'
l = partialModelSize;
noCol = size(matrix,2);
noRow = size(matrix,1);
count = 1;
for s=2:noCol
%s=3;
    % initialising the model for a value of 's', where 's' represents the number of cloumns present in the model
    model = [];
    for ta = 1:noCol
        for tb = 1:noCol
            if ta~=tb
                % Calculating the significant score for each element
                p = getScore(matrix,ta,tb,s);
                % Updating partial model
                model = updateModel(model,p,s,ta,tb,l+1,[],[],[],[],[],[]);
            end
        end
    end
    % if s==2, Realize and return best complete model of size 2 because the model is already complete and cannot be extended
    if s==2
        if size(model,2) == 0
            finModel = [];
        end
        bestIndex = 0;
        bestScore = -1;
        for i=2:size(model,2)
            if model(i).score > bestScore
                bestIndex = i;
                bestScore = model(i).score;
            end
        end
        [isInT,arrayTissue,prefinalmodel,isInG] = extend(matrix,[],model(bestIndex),[],[],[]);
        finModel = struct('GeneIndex',find(isInG),'TissueIndex',find(isInT));
        count = count + 1;                                          
    % else iterate until models are complete
    else                                                        
        InT=[];InG=[];
        for iter=0:s-3
            % Initialising the model
            extendedModel = struct('score',[],'noCol',[],'lowestRankCol',[],'highestRankCol',[],'noSmallElements',[],'noLargeElements',[],'extColumn',[],'modelNo',[]);   
            % Realizing best l models
            for i=1:size(model,2)-1                                   
                if iter == 0
                    [InT(i,:),arrayTissue(i,:),parent(i),InG(i,:)] = extend(matrix,[],model(i+1),[],[],[]);
                else
                    [InT(i,:),arrayTissue(i,:),parent(i),InG(i,:)] = extend(matrix,model(i+1).extColumn,model(i+1),arrayTissue(i,:),InT(i,:),InG(i,:));
                end
            end
            % Evaluating all possible extensions of the partial models until models are complete
            for i=1:size(parent,2)
                for t=1:noCol
                    if ~InT(i,t)
                        [extensionScore,ta,tb,arrayTissue(i,:)] = getScoreForExtension(matrix,t,arrayTissue(i,:),parent(i).noSmallElements,parent(i).noLargeElements,InG(i,:));
                        str = struct('score',extensionScore,'noCol',s,'lowestRankCol',ta,'highestRankCol',tb,'noSmallElements',parent(i).noSmallElements,'noLargeElements',parent(i).noLargeElements,'extColumn',t,'modelNo',i);  
                        [extendedModel,InG,InT] = updateModel(extendedModel,str.score,str.noCol,str.lowestRankCol,str.highestRankCol,l+1,parent(i).noSmallElements,parent(i).noLargeElements,t,InG,InT,i);
                    end
                end  
            end   
            clear model parent
            % Arranging the rows of 'InG', 'InT' and 'arrayTissue' matrices,
            % according to their respective order in the partial model
            mat1 = []; mat2 = []; mat3 = [];
            for iteration=2:size(extendedModel,2)
                Index = extendedModel(iteration).modelNo;
                mat1 = [mat1;InG(Index,:)]; mat2 = [mat2;InT(Index,:)];
                mat3 = [mat3;arrayTissue(Index,:)];
            end
            InG = mat1; 
            InT = mat2;
            arrayTissue = mat3;
            model = extendedModel;
        end
        % Finding the best complete model among all 'l' complete models
        bestIndex = 0;
        bestScore = -1;
        for i=2:size(extendedModel,2)
            if extendedModel(i).score > bestScore
                bestIndex = i;
                bestScore = extendedModel(i).score;
                isInT = InT(i-1,:);
                isInG = InG(i-1,:);
                bestarrayTissue = arrayTissue(i-1,:);
            end
        end
        newmodel = extendedModel(bestIndex);
        % Realizing the best model
        [isInT,arrayTissue,parent,isInG] = extend(matrix,newmodel.extColumn,newmodel,bestarrayTissue,isInT,isInG);
        % Checking for any overlap between newly added model and previous models
        if length(find(isInG))>1
            temp=[];
            for i=1:count-1
                if length(find(ismember(finModel(i).GeneIndex,find(isInG)))) == length(finModel(i).GeneIndex) && ...
                    length(find(ismember(finModel(i).TissueIndex,find(isInT)))) == length(finModel(i).TissueIndex)
                    temp=[temp,i];
                end
            end
            if isempty(temp)
                finModel(count) = struct('GeneIndex',find(isInG),'TissueIndex',find(isInT));
                count = count+1;
            else
                finModel(temp)=[];
                count = count-length(temp);
                finModel(count) = struct('GeneIndex',find(isInG),'TissueIndex',find(isInT));
                count = count+1;
            end
        end
        % Removing previous models which are smaller in size than newly
        % added model.
        %for i=1:count-2
          %  if size(finModel(i).GeneIndex,2)*size(finModel(i).TissueIndex,2) < size(finModel(count-1).GeneIndex,2)*size(finModel(count-1).TissueIndex,2)
            %    svModel = finModel(count-1);
              %  finModel(i:count-2)=[];
                %count=count-length(i:count-2);
                %finModel(count-1) = svModel;
                %break;
            %end
        %end
    end
clear model newmodel arrayTissue isInT isInG extendedModel parent ara InT InG;
end
end
