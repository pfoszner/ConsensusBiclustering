function result = infocontent(resfabic,matrix,nRow,nCol,p)
%% Function to calculate information content of biclusters
% 
% Inputs  :
%      resfabic  -  result of fabic.m (Structure consisting of factor matrix,loading 
%                   matrix, variational parameter matrix and noise variance
%        matrix  -  matrix obtained after preprocessing
%          nRow  -  number of rows in the matrix
%          nCol  -  number of coloums in the matrix
%             p  -  number of biclusters to be extracted
%
% Output  :
%        result  -  A structure consisting of the final loading matrix and factor matrix
% Author : Shruti Jain, 2014
%
% Contact : sjain@iitk.ac.in, sweetushruti963@gmail.com
%           Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India
%%

Z = resfabic.Z;
L = resfabic.L;
lapla = resfabic.lapla;
psi = resfabic.psi;
iin = 1/nCol;

vz = iin*(sum(Z.^2,2));
vz = sqrt(vz+10^-10);
ivz = 1./vz;

if size(ivz,1)==1
     nZ = ivz*Z;
     noL = vz*L;
     lapla = (vz^2)*lapla;
else
    nZ = bsxfun(@mtimes,Z',ivz')';
    noL = bsxfun(@mtimes,L,vz');
    lapla = bsxfun(@mtimes,lapla,(vz.^2)');
end

ini = zeros(nCol,p+1);
avini = zeros(1,p+1);
xavini = zeros(1,nCol+1);
idb = eye(p);
temp = bsxfun(@mtimes,noL',(1./psi)')';
ppL = noL'*temp;

for j = 1:nCol
    mat = idb+bsxfun(@rdivide,ppL',lapla(j,:))';
    ini(j,1:p) = log(diag(mat));
    s = log(det(mat));
    ini(j,p+1) = s;
    xavini(j) = s;
end

for i = 1:p
    avini(i) = sum(ini(:,i));
end

ss = sum(ini(:,p+1));
xavini(nCol+1) = ss;
avini(p+1) = ss;
Lz = noL*nZ;
M = eye(p);
U = matrix-Lz;

if avini(p+1)>10^-8 && p>1
    [~,index] = sort(avini(1:p),'descend');
    avini = avini(index);
    noL = noL(:,index);
    nZ = nZ(index,:);
    M = M(index,index);
end 

result =struct('noL',noL,'nZ',nZ);
end