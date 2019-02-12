function x = isafilter(x,t,dir)
% Matlab wrapper for call to C functions for itersa
%
% See Also: isaIterate, isaStep
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if strcmp(dir,'updown')
  misaBF_updown(x, t);
elseif strcmp(dir,'up')
  misaBF_up(x,t);
elseif strcmp(dir,'down')
  misaBF_down(x,t);
end
end