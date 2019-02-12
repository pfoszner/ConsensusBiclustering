function normalizedData = normalize(mat, type, error, maxit)
%Normalization
% Kluger Y., Basri R., Chang J.T. & Gerstein M., 
% Spectral Biclustering of Microarray Data: Coclustering Genes and Conditions
% Genome Research 2003.
%
% Inputs
%   matrix            - data matrix
%   normalization     - Normalization method to be applied to the matrix. As
%                       explained in Kluger et al, three methods are allowed:
%                       1)"log":Logarithmic Transformation
%                       2)"irrc":Independent Rescaling of Rows and Columns
%                       3)"bistochast":Bistochastization
%                     If "log" normalization is used, be sure you can apply
%                     logarithm to elements in data matrix, if there are 
%                     values under 1, it automatically will sum to each 
%                     element in mat (1+abs(min(mat))) Default is "log", 
%                     as recommended by Kluger et al.
%
% See Also: kSpectral
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if strcmp(type,'irrc')
  normalizedData=irrc(mat);
else
  if strcmp(type,'bistochast')
    normalizedData=bistochast(mat, error, maxit);
  elseif strcmp(type,'logt')
    normalizedData=logt(mat);
  end
end

end

%##########################################################################
function An = irrc(A)
% Independent Rescaling of Rows and Columns
% *Input*:
%    A     - matrix
R = sum(A,2).^(-0.5);
C = sum(A,1).^(-0.5);
An = ((bsxfun(@times, A, R))'.*C)'; % Error here
end
%##########################################################################

function An = bistochast(A, error, maxit)
%% Bistochastization 
% Basically a rescaling until convergence i.e. either maximum number of 
% iterations has been reached or change between two bistochastizations is 
% lesser than a threshold.
%
% *Input*:
%    A     - matrix to bistochastize
%    error - minimum change between two iterations of bistochastization.
%            Default:1e-9
%    maxit - maximum number of bistochastizations. Default: 10^3

if nargin<3
  error = 0.000000001;
  maxit = 1000;
end

diffr = error+1;
numit = 0;
An = A;

while (diffr>error && numit<maxit)
  Ap=An;
  R = sum(A,2).^(-0.5);
  C = sum(A,1).^(-0.5);
  An = (bsxfun(@times,R,Ap)'*C)'; % Error here
  diffr = abs(sum(An(:)-Ap(:)));
  numit = numit+1;
end
end
%##########################################################################

function K = logt(A)
%% Logarithmic Transformation and augmentation by additive inverse
% *Input*:
%    A - matrix to transform. Note: If A values are not greater than 1, log transform will be errorneous
n=size(A,1);
m=size(A,2);
L=log(A);
Li = sum(L,2)/m; %Sum rows
Lj = sum(L,1)/n; %Sum Columns
Lij = sum(L(:))/(n*m);
a1 = bsxfun(@minus,L,Li); %L-Li
a2 = (bsxfun(@plus,a1,Lij));%(L-Li+Lij)'
K = (bsxfun(@minus,a2,Lj)); % ((L-Li+Lij)'-Lj)'
end
