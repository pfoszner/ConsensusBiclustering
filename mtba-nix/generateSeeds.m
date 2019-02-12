function g = generateSeeds(len, count, method, sparcity)
%Generate random input seeds for ISA
%
% Usage 
% >> g = generateSeeds(length, count, method, sparcity)
%
% Inputs
%   length    - the length of the seeds, should be the number of rows in
%               your input data for row seeds and number of columns for 
%               column seeds.
%   count     - the number of seeds to generate
%   method    - the method for generating the seeds. Right now only 'uni'
%               is supported. It picks one element each seed uniformly but
%               randomly.
%   sparcity  - integer number giving the number of non-zero values in each
%               seed vector which is then recycled to have tge same length 
%               as the number of seeds.
%
% Output
%   g         - numeric matrix with 0/1 value
%
% See Also: itersa, isafilter, isaFilterRobust, isaIterate, isaNormalize
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<3
  count     = 100;
  method    = 'uni';
  sparcity  = 2;
end

if strcmp(method, 'uni')
  sparcity = repmat(sparcity, count,1);
  g = zeros(len, count);
  for i=1:count
    g(randsample(len,sparcity(i)),i) = 1;
  end
end
end