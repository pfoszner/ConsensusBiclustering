function biClustResult = fuzzyBiclust(A)
%matlab code : A modified Fuzzy Co-clustering (MFCC) Approach for
%                microarray data analysis
%based on : A partitioning based algorithm to fuzzy co-cluster documents
%              and words (FCR)
%Paramaters to be given as input :
%  A=input matrix
%  C=no of biclusters to be formed
%  tu,tv = pre-defined degree of fuzziness parameters
%  epsilon = parameter for convergence
%  alpha = threshold parameter for fuzzy membership
% Author - Esha Dutta, Summer Intern, Intelligent Informatics lab, IITK
clearvars -except A;
N = size(A,1);  %number of genes
K = size(A,2);  %number of conditions
C = input('enter the number of biclusters to be formed :');
disp(sprintf('DEFAULT VALUE OF:\nepsilon = 0.01\nthreshold = 0.5\npre-defined degree of fuzziness for genes = 0.7\npre-defined degree of fuzziness for conditions = 0.7 '));
ch=input('want to change the values?(yes=1/no=0):');
epsilon =0.01;
alpha =0.5;
tu =0.7;
tv =0.7;
if ch == 1
epsilon = input('enter epsilon :');
alpha = input('enter threshold :');
tu = input('pre-defined degree of fuzziness for genes :');
tv = input('pre-defined degree of fuziness for conditions :');
else
for i=1:N
    a(:,i)=rand(C,1);
    b(i)=sum(a(:,i));
    if b(i) == 0
        disp('incorrect initial membership generated, run the code once again to obtain output');
        flag = false;
    else
    u(:,i)=a(:,i)/b(i);
    end
end
flag=true;
while (flag)
    v = fuzzyBiclust_compute_v(u,A,N,K,C,tv);
    u1 = fuzzyBiclust_compute_u(v,A,N,K,C,tu);
    if fuzzy_check1(u,u1,epsilon)==1
        flag=false;
    else
        u=u1;
        flag=true;
    end
    end
u=fuzzy_modified(u1,alpha);
v=fuzzy_modified(v,alpha);
v=v';
biClustResult.RowxNum = logical(u);
biClustResult.NumxCol = logical(v);
biClustResult.ClusterNo = C;
for k = 1:C
    rows = find(biClustResult.RowxNum(:,k)>0);
    cols =find(biClustResult.NumxCol(k,:)>0);
    biClustResult.Clust(k) = struct('rows',rows,'cols',cols);
end
end