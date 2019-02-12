function rez=LAS_score(sum, subrows, subcols, mrows, mcols)
%Function to compute the score of potential biclusters
% The computation of score is based on Bonferroni significance correction. 
% The score function trades off between the submatrix size and average value.
% Input:
%        sum      :   Cumulative sum of potential biclusters 
%        subrows  :   Subset of rows
%        subcols  :   Subset of columns
%        mrows    :   Number of rows in input matrix
%        mcols    :   Number of column in input matrix
% Output:
%        rez      :   Calculated score for potential bicluster 
% See Also: LAS
% -------------------------------------------------------------------------
% Modified version of LAS_score.m in LAS_Matlab package by Andrey Shabalin
% Changes required to include the program in the toolbox.
% 2013 Sumanik Singh <sumaniksingh@gmail.com>
%--------------------------------------------------------------------------

cnrows = gammaln(mrows+1) - gammaln(subrows+1) - gammaln(mrows-subrows +1);
cncols = gammaln(mcols+1) - gammaln(subcols+1) - gammaln(mcols-subcols +1);
ar = sum./sqrt(subrows.*subcols);
rest2 = - ar.^2/2 + log( erfcx(ar/sqrt(2))/2 );
rez = - rest2 - cnrows - cncols;
end
