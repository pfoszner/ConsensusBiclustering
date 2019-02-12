function [p,ta,tb,arrayTissue] = getScoreForExtension(dataset,extColumn,arrayTissue,noSmallElements,noLargeElements,isInG)
% Function to compute significance score for extended model
% The function is modified version of 'getScore' function, which computes
% for significance score for partial models.
% Inputs:
%        dataset          :    Input matrix
%        extColumn        :    Column being considered for extending the partial model
%        arrayTissue      :    Array containing column numbers present in the models 
% If partial model is of order (a,b), then 
%        noSmallElements  :    a, i.e number of smaller elements 
%        noLargeElements  :    b, i.e number of larger elements 
%        isInG            :    (i,j) entry is 1, if i-th model contains j-th row 
% Outputs:
%        p                :    Significance score for the extended model
%        ta               :    Column with lowest rank
%        tb               :    Column with highest rank
%        arrayTissue      :    Array containing column numbers present in the models
%
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of
%          Technology, Kanpur, India
%% Main routine for the function
% Assigning values to the parameters required for implementation
a = noSmallElements;
b = noLargeElements;
t = extColumn;
noRow = size(dataset,1);
noCol = size(dataset,2);
% Finding ranked matrix D
[so,ind] = sort(dataset,2);
[s1,D] = sort(ind,2);
ta = arrayTissue(a);
tb = arrayTissue(length(arrayTissue)-b+1);
% The if condition calculates exact value for p if the extended model is already complete
if a+b==length(arrayTissue)-1
    k=0;
    % Getting number of compatible rows
    for row=1:noRow
        if (isInG(row) && D(row,ta)<D(row,t) && D(row,t)<D(row,tb))
            k=k+1;
        end
    end
    p=k/noRow;
% else condition is used when the model is still not complete. The
% significance score is estimated using newton iteration method to solve equation.
else
    commonDenominator = nchoosek(noCol,length(arrayTissue));
    % Computing B[i] or B
    B = 1/(nchoosek(noCol,a+b+1)*factorial(a+b+1));
    % Computing A[i] and diff[i] = A[i]-B
    for i=1:noRow
        A(i) = 0;
        if (isInG(i) && D(i,ta)<D(i,t) && D(i,t)<D(i,tb))
            gi = 0;
            if a==b
                gi = D(i,tb) - D(i,t) - 1;
            else
                gi = D(i,t) - D(i,ta) - 1;
            end
            if gi >= length(arrayTissue)-a-b-1 && gi>0 && (length(arrayTissue)-a-b-1)>0
                A(i) = nchoosek(gi,length(arrayTissue)-a-b-1)/commonDenominator;
            else
                A(i)=0;
            end
        end
            diff(i) = A(i) - B;
    end
    % Starting Point for Newton Iteration
    p = 0.05;              
    % Performing 20 newton iteration steps to estimate p 
    for iteration = 1:20
         f = -noRow;
         dfdp = 0;
         for i = 1:noRow
             denominator(i) = diff(i)*p + B;
             if denominator(i)==0 && A(i)~=0
                 disp('WARNING: Division by 0 in Newton Iteration.');
                 f = 0;
                 dfdp = 1;
             else if A(i)~=0
                 f = f + (A(i)/denominator(i));
                 dfdp = dfdp - (A(i)*diff(i)/(denominator(i)^2));
                 end
             end
         end
         if dfdp==0
             p = 0;
         else
             q = p;
             p = p - (f/dfdp);
             if abs(f/dfdp) < 0.0001
                 break;
             end
             if p<0 && iteration<20
                 p = q/2;
             end
         end
    end
    % special modification of newton algorithm, required to prevent from being caught in negative range 
    if p<0 || p>1
         p=0;
    end
end 
end
