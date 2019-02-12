function [res] = fuzzy_check1( x,y,e )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
pe=abs(x-y);
d=max(max(pe));
if d<=e
  res=1;
else
    res=0;
end
    
end

