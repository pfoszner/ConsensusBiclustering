function result = compare_row(row1,row2,coloum)
% Function gives result as coloum vector of the indices for which nonzero
% integer values of both the row is same
%
% Inputs :
%   row1   : row1 to be compared
%   row2   : row2 to be compared
%   coloum : coloum vector along which we have to compare
%
% Outputs :
%   result : row vector of the indices
%
%% Main program
 n = size(row1,2);
 result = struct('res',[]);
 result.res = zeros(1,n);
for i = 1:n
    if coloum(i) == i
        if row1(i) == row2(i) &&row1(i)~=0
            result.res(i) = 1;
        end
    end
end

end
        


