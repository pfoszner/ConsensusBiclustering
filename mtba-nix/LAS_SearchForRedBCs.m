function TheBiclusters = LAS_SearchForRedBCs(data, numBCfind, iterationsPerBC,scoreThreshold,option)
%% Function to search for biclusters in input matrix
% This function performs all the comparisons and calculations relating to
% finding biclusters.
% Input:
%        data             :    Input matrix whose biclusters are to be found
%        numBCfind        :    Number of biclusters to be found
%        iterationsPerBC  :    Number of iterations per Bi-cluster
%        scoreThreshold   :    Threshold value of score. If the value of score
%                              s below the threshold value then the iterations are 
%                              stopped.
%        option           :    String containing information about data
%                              (either 'positive' or 'negative')
% Output:
%        TheBiclusters    :    This variable stores the best bicluster found
% -------------------------------------------------------------------------
% Modified version of LAS_SearchForRedBCs.m in LAS_Matlab package by Andrey Shabalin
% Changes required to include the program in the toolbox.
% 2013 Sumanik Singh <sumaniksingh@gmail.com>
%--------------------------------------------------------------------------

	TheBiclusters = struct('score', {}, ...
				'rows', {}, ...
				'cols', {}, ...
				'avg', {});
    disp(['Biclusters obtained for ' option ' data']);
	for bcNum = 1:numBCfind
		bestBC.score = -Inf;
		parfor i=1:iterationsPerBC
			bc = SingleBiclusterSearch(data);
			bestBC = BetterBicluster(bestBC,bc);
		end;
		disp(['BC#=' num2str(bcNum) ' found,	' num2str(sum(bestBC.rows)) 'x' num2str(sum(bestBC.cols)) ',	score=' num2str(bestBC.score)]);
		TheBiclusters(bcNum) = bestBC;

		data(bestBC.rows,bestBC.cols) = ...
			data(bestBC.rows,bestBC.cols) - bestBC.avg;

		if(bestBC.score < scoreThreshold)
			break;
		end;
	end;
end

function bc = BetterBicluster(bc1, bc2)
%% Function to compare the scores of two biclusters
% The function compares the scores of two biclusters found and returns the
% better among them.
% Input:
%        bc1   :  Bicluster 1
%        bc2   :  Bicluster 2
% Output:
%        bc    :  Better bicluster
%% Comparing the scores
	if (bc1.score>bc2.score)
		bc = bc1;
	else
		bc = bc2;
	end;
end
