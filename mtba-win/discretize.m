function dismat = discretize(x, nof, quant)
%Discretize the matrix values into bins
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<2
  nof =10;
  quant = false;
end

dismat = x;
[ni nj] = size(x);

if quant
  levels = quantile(x, nof-1);
else
  mindata = min(x(:));
  maxdata = max(x(:));
  levels = zeros(nof,1);
  diff = (maxdata-mindata)/nof;
  
  for k=1:nof
    levels(k) = mindata + k*diff;
  end
end

for i=1:ni
  for j=1:nj
    for k=1:nof
      if x(i,j)<=levels(nof-k+1), dismat(i,j) = k; end
    end
  end
end
