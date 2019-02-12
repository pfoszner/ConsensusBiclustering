function [rowSets, colSets, numBiclusters] = LASmainfile(data, BCnumPositive, BCnumNegative, iterationsPerBC, scoreThreshold)
%[rowSets, colSets, numRedBiclusters, numGrnBiclusters] = LASmainfile(data,
%BCnumPositive, BCnumNegative, iterationsPerBC, scoreThreshold)
%% Implementation of LAS
% Inputs:
%        data           :     Input matrix whose biclusters are to be found
%        BCnumPositive  :     Number of Bi-clusters for positive data i.e 
%                (data considered without changing the sign, InputMatrix = data)
%        BCnumNegative  :     Number of Bi-clusters for negative data i.e 
%              (data considered after reversing the sign, InputMatrix = -(data))
%        iterationPerBC :     Number of iterations per Bi-cluster
%        scoreThreshold :     Threshold value of score. If the value of score
%                             s below the threshold value then the iterations are 
%                             stopped.
% Outputs:
%        rowSets          :     It is a logical matrix containining information about rows
%                               contained in different biclusters.
%                               If (i,j) value is 1, then j row is contained in i bicluster.
%        colSets          :     It is a logical matrix containining information about columns
%                               contained in different biclusters.
%                               If (i,j) value is 1, then j column is contained in i bicluster.
%        numRedBiclusters :     Number of biclusters for positive data.
%        numGrnBiclusters :     Number of biclusters for negative data.
%
% -------------------------------------------------------------------------
% Modified version of main.m in LAS_Matlab package by Andrey Shabalin
% Changes required to include the program in the toolbox.
% 2013 Sumanik Singh <sumaniksingh@gmail.com>
%--------------------------------------------------------------------------

%% preprocessing
% Column Center
CC = @(x)(x-repmat(mean(x,1),[size(x,1) 1]));

% Column Standardize
CS_ = @(x)(x./repmat(std(x),size(x,1),1));
CS = @(x)CS_(CC(x));

%  The bend tails LN transform
LN_ = @(x)(sign(x).*log(1+abs(x)));
LN = @(x)CS(LN_(x));

data = CS(data);
dataLN = LN(data);

dataCur = mean(data(:).^4);
dataLNcur = mean(dataLN(:).^4);

if(abs(dataLNcur-3) < abs(dataCur-3))
	data = dataLN;
end;

clear dataLN;
clear RMC CC CS CS_ LN LN_ dataCur dataLNcur;
%% Calling the main program
RedBiclusters = LAS_SearchForRedBCs( data,BCnumPositive,iterationsPerBC,scoreThreshold,'positive');
GrnBiclusters = LAS_SearchForRedBCs(-data,BCnumNegative,iterationsPerBC,scoreThreshold,'negative');

%% Saving the results
arrayScores = [];
for i=1:length(RedBiclusters)
    arrayScores = [arrayScores,RedBiclusters(i).score];
end
for i=1:length(GrnBiclusters)
    arrayScores = [arrayScores,GrnBiclusters(i).score];
end

[score, index] = sort(arrayScores,'descend');
rowSetsCell = cell(BCnumPositive,1);
colSetsCell = cell(BCnumPositive,1);

if BCnumPositive < (length(RedBiclusters) + length(GrnBiclusters))
    numBiclusters = BCnumPositive;
    for i=1:numBiclusters
        if(index(i)>length(RedBiclusters))
            rowSetsCell{i} = GrnBiclusters(index(i)-length(RedBiclusters)).rows';
            colSetsCell{i} = GrnBiclusters(index(i)-length(RedBiclusters)).cols;
        else
            rowSetsCell{i} = RedBiclusters(index(i)).rows';
            colSetsCell{i} = RedBiclusters(index(i)).cols;
        end
    end
else
    numBiclusters = length(RedBiclusters) + length(GrnBiclusters);
    for i = 1:numBiclusters
        if(index(i)>length(RedBiclusters))
            rowSetsCell{i} = GrnBiclusters(index(i)-length(RedBiclusters)).rows';
            colSetsCell{i} = GrnBiclusters(index(i)-length(RedBiclusters)).cols;
        else
            rowSetsCell{i} = RedBiclusters(index(i)).rows';
            colSetsCell{i} = RedBiclusters(index(i)).cols;
        end
    end
end

rowSets = cat(1,rowSetsCell{:});
colSets = cat(1,colSetsCell{:});
end
