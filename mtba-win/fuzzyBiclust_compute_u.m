function y = fuzzyBiclust_compute_u(p,q,n,k,c,t1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a=n*k;
sum1=0;  %sum of all elements of q
sum1=sum(sum(q));
p1=sum1/a;
a1=p*q';
for i=1:c
    s1(i)= sum(p(i,:));
end
for i=1:c
    for j=1:n
        v1(i,j)=a1(i,j)-(p1*s1(i));
        v2(i,j)=v1(i,j)/t1;
        v3(i,j)=exp(v2(i,j));
    end
end
for i=1:c
   for j=1:n
        y(i,j)=v3(i,j)/sum(v3(:,j));
   end
end
    

end

