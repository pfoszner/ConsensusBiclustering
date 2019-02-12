function matd = binarize(datamatrix, threshold)
%Methods to convert a real matrix to a binary matrix as per median 
%
% Usage
% >> binarizedmatrix = binarize(datamatrix, threshold)
%
% Inputs: 
%   datamatrix            - The data matrix to be binarized
%   threshold (optional)  - Threshold used to binarize. 
%                           Values over threshold will be set to 1, the rest to 0.
%   
% Outputs:
%   matd    - Binarized matrix
%   
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin < 1
    error('input :  No matrix as input');
end

matd = datamatrix;

if nargin < 2
    %# Calculate Median of entire matrix as threshold
    threshold = min(datamatrix(:)) + (max(datamatrix(:)) - min(datamatrix(:)))/2;
    fprintf('Threshold: %d\n', threshold);
end

matd(matd<=threshold) = 0;
matd(matd>threshold)  = 1;

end