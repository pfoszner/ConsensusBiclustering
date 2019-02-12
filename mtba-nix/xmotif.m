function biClustResult = xmotif(matrix, ns, nd, sd, alpha, number)
%XMotif Biclustering
% Extracting conserved gene expression motifs from gene expression data
% Murali TM, Kasif S
% Bioinformatics Program, 48 Cummington St., Boston University, Boston, MA 02152, USA.
% Pac Symp Biocomput. 2003;:77-88.
%
% Usage
% >> biClustResult = xmotif(matrix, ns, nd, sd, alpha, number)
%
% Inputs:
%   matrix            - the data matrix (discrete)
%   ns                - Number of rows to be chosen
%   nd                - Number of repetitions
%   sd                - Sample size in repetitions
%   alpha             - Scaling factor for column result
%   number            - Number of bicluster to be found
%
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] if Row i is 
%                     in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] if Col j is 
%                     in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices.
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

%##########################################################################
%       Error Checking and Default Values
%##########################################################################
if nargin < 1
  error('input :  No matrix as input');
end

if nargin < 2
  ns      = 10;
  nd      = 10;
  sd      = 5;
  alpha   = 0.05;
  number  = 100;
end
%##########################################################################

mat   =  matrix;
nrow  = size(matrix, 1);
ncol  = size(matrix, 2);
x     = false(nrow, number);
y     = false(number, ncol);

STOP  = false;
logr  = true(nrow,1);

for i = 1:number
  if sum(logr)<2
    STOP = true;
    break;
  end
  
  erg = bigxmotif(mat,ns,nd,sd,alpha);
  
  if sum(erg.gen)==0
    STOP = true;
    break;
  else
    x(logr,i)   = erg.gen;
    y(i,:)      = erg.co;
    logr(logr)  = ~(logr(logr)&erg.gen);
    mat         = matrix(logr,:);
  end
end

if STOP
  biClustResult.RowxNum   = x(:,1:(i-1));
  biClustResult.NumxCol   = y(1:(i-1),:);
  biClustResult.ClusterNo = i-1;
  for j=1:(i-1)
    rows = find(biClustResult.RowxNum(:,j)>0);
    cols =find(biClustResult.NumxCol(j,:)>0);
    biClustResult.Clust(j) = struct('rows', rows, 'cols', cols); %
  end
else
  biClustResult.RowxNum   = x;
  biClustResult.NumxCol   = y;
  biClustResult.ClusterNo = i;
  for j=1:i
    rows = find(biClustResult.RowxNum(:,j)>0);
    cols =find(biClustResult.NumxCol(j,:)>0);
    biClustResult.Clust(j) = struct('rows', rows, 'cols', cols); %
  end
end
end

%##########################################################################
%   FindMotif() algorithm from the paper
%##########################################################################
function res = bigxmotif(mat,ns,nd,sd,alpha)
% See something about preprocessing

siz   = 4;
nc    = size(mat,2); % Column Size
gen   = false(size(mat,1),1);
co    = false(size(mat,2),1);

for i=1:ns
  ci        = randsample(1:nc,1);
  logc      = true(size(mat,2),1);
  logc(ci)  = false;
  
  for j=1:nd
    D     = randsample(1:nc,sd,true,logc);
    gci   = mat(:,ci);
    gciD  = [D ci];
    rS    = sum(bsxfun(@eq,mat(:,gciD),gci),2); % Row Sum
    gij   = bsxfun(@eq,rS,size(gciD,2));
    
    if (sum(gij)>=max([sum(gen) 2]))
      cci = mat(gij,ci);
      cS  = sum(bsxfun(@eq,mat(gij,:),cci),1); % Column sum
      cij = bsxfun(@eq,cS,sum(gij));
      
      if (sum(cij)>=(alpha*nc) && ((sum(gij)*sum(cij))>siz))
        gen = gij;
        co = cij;
        siz = sum(gij)*sum(cij);
      end
    end
  end
end
res = struct('gen',gen,'co',co);
end

% function res = storexmotif(mat,ns,nd,sd,alpha)
% nc = size(mat,2);
% xstore = zeros(size(mat,1),ns*nd);
% ystore = zeros(ns*nd,size(mat,2));
% for i=1:ns
%     ci = randsample(1:nc,1);
%     logc = true(size(mat,2),1);
%     logc(ci)=false;
%
%     for j=1:nd
%         D = randsample(1:nc,sd,true,logc);
%
%         gci = mat(:,ci);
%         gciD = [D ci];
%         rS = sum(mat(mat(:,gciD)==gci),2); % Row Sum
%         gij = bsxfun(@eq,rS,size(gciD,2));
%
%         if sum(gij)>2
%             cci = mat(gij,ci);
%             cS = sum(mat(mat(gij,:)==cci),1); % Column sum
%             cij = bsxfun(@eq,cS,sum(gij));
%
%             xstore(:,j+(i-1)*nd)=gij*1;
%             ystore(j+(i-1)*nd,:)=cij*1;
%         end
%     end
% end
% res = struct('xstore',xstore,'ystore',ystore);
% end
