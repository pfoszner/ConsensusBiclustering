function bc = SingleBiclusterSearch(data)
%% Search for a bicluster
% This function searches through dataset for biclusters by firstly searching 
% with fixed 'bcNumRows' and'bcNumCols' and then with both variables being flexible.
% Input:
%        data    :   Input matrix whose biclusters are to be found
% Output:
%        bc      :   Bicluster found
% -------------------------------------------------------------------------
% Modified version of SingleBiclusterSearch.m in LAS_Matlab package by Andrey Shabalin
% Changes required to include the program in the toolbox.
% 2013 Sumanik Singh <sumaniksingh@gmail.com>
%--------------------------------------------------------------------------

[m n] = size(data);
bcNumRows = floor(rand(1)^2*m/2)+1;
bcNumCols = floor(rand(1)^2*n/2)+1;

% start with a random columnset
randomColPermutation = randperm(n);
colSet = (randomColPermutation<=bcNumCols);
clear fraction2select randomColPermutation;
%% First we search with bcNumRows and bcNumCols fixed
prevAvg = -Inf;
currAvg = 0;
while(prevAvg~=currAvg)
	prevAvg = currAvg;
	
	% calculate the row sums over the selected columns
	rowSums = data*colSet';

	% sort sums, saving the info about order
	[sortedRowSums orderRowSums] = sort(rowSums,'descend');

	% select bcNumRows rows with larges averages
	rowSet = zeros(m,1);
	rowSet(orderRowSums(1:bcNumRows)) = true;
	
	% calculate the column sums over the selected rows
	colSums = rowSet'*data;

	% sort sums, saving the information about permutation]
	[sortedColSums orderColSums] = sort(colSums,'descend');

	% select bcNumRows rows with larges averages
	colSet = zeros(1,n);
	colSet(orderColSums(1:bcNumCols)) = true;
	currAvg = mean(sortedColSums(1:bcNumCols)/bcNumRows);
end;
clear prevAvg currAvg rowSums colSums sortedRowSums sortedColSums;
clear orderRowSums orderColSums;
%% now let bcNumRows and bcNumCols be flexible
prevScore = -Inf;
currScore = 0;
while(prevScore~=currScore)
	prevScore = currScore;
	
	% calculate the row sums over the selected columns
	rowSums = data*colSet';

	% sort sums, saving the info about order
	[sortedRowSums orderRowSums] = sort(rowSums,'descend');

	% Calculate the sums of potential biclusters
	potBcSumsR = cumsum(sortedRowSums);
    
	% the Scores of the potential biclusters
	potScoresR = LAS_score(potBcSumsR,(1:m)',bcNumCols,m,n);
	[maxPotScoreR bcNumRows] = max(potScoresR);	
	rowSet = zeros(m,1);
	rowSet(orderRowSums(1:bcNumRows)) = true;

	% calculate the column sums over the selected rows
	colSums = rowSet'*data;

	% sort sums, saving the information about permutation]
	[sortedColSums orderColSums] = sort(colSums,'descend');

	% Calculate the sums of potential biclusters
	potBcSumsC = cumsum(sortedColSums);
	
    % the Scores of the potential biclusters
	potScoresC = LAS_score(potBcSumsC,bcNumRows,(1:n),m,n);
	[maxPotScoreC bcNumCols] = max(potScoresC);	
	
    % select bcNumRows rows with larges averages
	colSet = zeros(1,n);
	colSet(orderColSums(1:bcNumCols)) = true;
	currScore = maxPotScoreC;
end;
bc = struct( ...
	'score', currScore, ...
	'rows', logical(rowSet), ...
	'cols', logical(colSet), ...
	'avg', mean(sortedColSums(1:bcNumCols))/bcNumRows);
return;
