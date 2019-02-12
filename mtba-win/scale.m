function scaled = scale(data)
%Scale the matrix to get z scores column wise
%
% See Also: isaIterate, isaNormalize
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

% First mean center the columns
centered = detrend(data, 'constant');
% Divide by std. deviation
sigma = std(centered);
scaled = bsxfun(@rdivide, centered, sigma);
end
