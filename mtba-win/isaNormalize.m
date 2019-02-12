function normed = isaNormalize(data, prenormalize)
%Normalize the data matrix getting their z scores
%
% Usage
% >> normed = isaNormalize(data, prenormalize)
%
% Inputs
%   data         - the input numeric matrix
%   prenormalize - true/false
%                   If true, then row-wise scaling is calculated on the 
%                   columnwise scaled matrix and not on the input matrix 
%                   directly
%
% Output
%   normed      - normalized data consisting of two matrices, first one
%                 transposed. Also contains information if the data was
%                 prenormalized or contains NaN.
%
% See Also: itersa, isaFilterRobust, scale
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<1
  error('input : No matrix as input')
end

if nargin<2
  prenormalize = false;
end

if ~ismatrix(data)
  error('input : `data` must be a matrix');
end

% If OK, Normalize
Ec = scale(data');

if prenormalize
  Er = scale(Ec');
else
  Er = scale(data);
end

normed = struct('Er',Er','Ec',Ec','prenormalize',prenormalize,'hasNaN',any(isnan(Er(:)))||any(isnan(Ec(:))));
end