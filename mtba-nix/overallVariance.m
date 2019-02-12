function ov = overallVariance(x)
% Gives overall variance as mean of row and column variance
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

rv = variance(x,2); % Row
cv = variance(x,1); % Col
[n m] = size(x);

ov = (n.*rv+m.*cv)./(n+m);
end