function result = preprocessing(matrix,p,alpha,cyc,spl,spz,non_negative,random,center,norm,scale,lap,nL,lL,bL)
%% Function for preprocessing and initialization of all variables.
% Also calls fabic.m for further processing.
%
% Inputs :
%             matrix  -  Input datamatrix
%                  p  -  Number of Biclusters to be found;default = 13
%              alpha  -  sparseness loadings;default = 0.01
%                cyc  -  Number of Iterations;default = 500
%                spl  -  sparseness prior loadings;default = 0;
%                spz  -  sparseness factor;default = 0.5
%       non_negative  -  non_negative loading and factors :0 for no, 1 for yes:default = 0
%             random  -  initialization of loadings : <=0 for svd, >0 for random initialization;default = 1
%             center  -  data centering : 1 for mean, 2 for median, 3 for mode and 0 for no centering;default = 2
%               norm  -  data normalization : 1 for (0.75-0.25 quantile), >1 for (var=1), 0 for no normalization;default = 1
%              scale  -  loading vectors are scaled in each iteration to the given variance. 0.0 indicates non scaling; default = 0        
%                lap  -  minimal value of the variational parameter, default = 1
%                 nL  -  maximal number of biclusters at which a row element can participate; default = 0(no limit)
%                 lL  -  maximal number of row elements per bicluster; default = 0 (no limit)
%                 bL  -  cycle at which the nL or lL maximum starts; default = 0 (start at the beginning)
%  
% See Also : fabic.m
% 
% Output : 
%             result  -  Final Bicluster Model
%
% Author: Shruti jain, 2014
%        
% Contact: sjain@iitk.ac.in, sweetushruti963@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India
%% preprocessing the data

nRow = size(matrix,1);
nCol = size(matrix,2);
if p > min(nCol,nRow)
    error('Too many biclusters.');
end

init_lapla = 1;
init_psi = 0.2;
iin = 1/nCol;

cent = zeros(nRow,1);

if center == 1
    cent = mean(matrix,2);
end
if center == 2
    cent = median(matrix,2);
end
if center == 3
    cent = estimateMode(matrix);
end

matrix = bsxfun(@minus,matrix,cent);
xx = ones(nRow,1);
scaleData = ones(nRow,1);

if norm > 0
    if norm > 1
        scaleData = 1./sqrt(iin*(sum(matrix.*matrix,2))+0.001*ones(1,nRow)');
    else
        for i = 1:nRow
            scaleData(i) = quant(sort(matrix(i,:)),nCol,0.75)-quant(sort(matrix(i,:)),nCol,0.25);
            scaleData(i) = sqrt(scaleData(i)+0.001);
        end
        scaleData = 1./scaleData;
    end
    matrix = bsxfun(@mtimes,matrix',scaleData')';
end

if random > 0
    L = randn(nRow,p);
else
    [u,s,~] = svd(matrix);
    L = u(:,1:p)*s(1:p,1:p);
end

if scale > 0
    dL = 1./sqrt(sum(L.*L,2)+0.001*ones(1,nRow)'); 
    L = bsxfun(@mtimes,L',(scale*dL)')';
end

lapla = init_lapla*ones(nCol,p);
psi = init_psi*xx;

result = fabic(matrix,psi,L,lapla,cyc,alpha,spl,spz,scale,lap,nL,lL,bL,non_negative);

end


%% Function to evaluate quantile of a row vector  %
function result = quant(row,nCol,f)

index = f*(nCol-1);
lhs = floor(index);
delta = index-lhs;
if delta ==0
   result = (1-delta)*row(lhs+1);
else
    result = (1-delta)*row(lhs+1) + delta*row(lhs+2);
end
end


%% Function to estimate mode of the matrix
function result = estimateMode(matrix,maxiter,tol,alpha,a1,G1)

maxiter = 50;
tol = 0.001;
alpha = 0.1;
a1 = 4;
G1 = 0;

nRow = size(matrix,1);
nCol = size(matrix,2);
iin = 1/nCol;
u = median(matrix,2);
xu = bsxfun(@minus,matrix,u);
xx = sqrt(iin*(sum(xu.*xu,2))+0.001*ones(1,nRow)');
dxx = 1./xx;
matrix = bsxfun(@mtimes,matrix',dxx')';

u = median(matrix,2);
gxu = 10*tol;
iter = 1;
alpha = alpha/nCol;
xu = bsxfun(@minus,matrix,u);

if G1
    while abs(gxu(1))>tol && iter<maxiter
        gxu = alpha*(sum(tanh(a1*xu),2));
        u = u+gxu;
        xu = bsxfun(@minus,matrix,u);
        iter = iter+1;
    end
else
    while abs(gxu(1))>tol && iter<maxiter
        gxu = alpha*(sum(xu.*exp(-a1*(xu.*xu)/2),2));
        u = u+gxu;
        xu = bsxfun(@minus,matrix,u);
        iter = iter+1;
    end
end
u = xx.*u;
result = u;
end
