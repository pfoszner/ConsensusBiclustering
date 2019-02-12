function y = fuzzy_modified( a,b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y=a';
p=size(y,1);
q=size(y,2);
for i=1:p
    for j=1:q
        if y(i,j)>=b
            y(i,j)=1;
        else
            y(i,j)=0;
        end
    end
end

            

end

