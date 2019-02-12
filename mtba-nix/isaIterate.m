function isaresult = isaIterate(normed_data, row_seeds, col_seeds, ...
                                thr_row, thr_col, direction, maxiter,...
                                convergence, cor_limit, eps, corx)
%Perform ISA on the (normalized) input matrix.
%
% Usage
% >> isaresult = isaIterate(normed_data, row_seeds, thr_row, ...
%                               col_seeds, thr_col, direction,...
%                               convergence, cor_limit, eps, corx, maxiter)
%
% Inputs
%   normed_data     - the normalized data consisting of two matrices, first
%                     one transposed. Also contains information if the data
%                     was prenormalized or contains NaN.
%   row_seeds       - row seed vectors to start
%   col_seeds       - column seed vectors to start
%   thr_row         - Threshold parameter/s for the rows. Higher values 
%                     indicate a more stringent threshold and the result 
%                     biclusters will contain less rows on average. The
%                     threshold is measured by the number of standard 
%                     deviations from the mean, over the values of the row 
%                     vector. If it is a vector then it must contain an 
%                     entry for each seed.
%   thr_col         - Threshold parameter/s for the columns. Analogue of
%                     row_seeds for columns
%   direction       - String cell array of length two, one for rows, second 
%                     for columns. Specifies whether we are interested in 
%                     rows/columns that higher('up') than average, lower 
%                     ('down') than average, or both ('updown').
%   maxiter         - maximum number of iterations allowed.
%   convergence     - Convergence criteria for ISA iteration.
%                    'cor':  convergence is measured based on high Pearson 
%                            correlation of the subsequent rows and columns
%                    'eps':  n all entries of the subsequent row and column
%                            vectors must be close to each other. Also see
%                            `eps` argument below.
%                    'corx': similar to 'cor', but the current row/column 
%                             vectors are compared to the ones corx steps 
%                             ago, instead of the ones in the previous
%                             step. Also see `corx` argument below.
%   cor_limit       - the correlation limit for convergence if 'cor' method
%                     is used.
%   eps             - limit of convergance if 'eps' method is used.
%   corx            - the number of iterations to use as a shift, for
%                     checking convergence with 'corx' method.
%
%
% Output
%   isaresult: A structure consisting of
%     Rows          - the row components of biclusters, every column
%                     corresponding to a bicluster.
%     Columns       - the column components of biclusters, in the same
%                     format as rows.
%     rundata       - structure with information about ISA runs.
%     seeddata      - matrix with information about the biclusters.
%
% See Also: isaFilterRobust, isaStep, isafilter, scale
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India


if isempty(row_seeds)&&isempty(col_seeds)
  error('No seeds, nothing to do');
end

if ~isempty(thr_row)&&isempty(thr_col)
  thr_col = thr_row;
end

if nargin==8
  if strcmp(convergence, 'eps')
    eps = 0.0001;
    corx = NaN;
    cor_limit = NaN;
  elseif strcmp(convergence, 'cor')
    cor_limit = 0.99;
    eps = NaN;
    corx = NaN;
  elseif strcmp(convergence, 'corx')
    corx = 3;
    eps = NaN;
    cor_limit = 0.99;
  end
end

if nargin<8
  convergence = 'corx';
  corx = 3;
  eps = NaN;
  cor_limit = 0.99;
end

if nargin<7
  maxiter = 100;
end

if nargin<6
  direction = {'updown', 'updown'};
end
  

if size(row_seeds,1)~=size(normed_data.Er,2)
  error('Invalid row seed length');
end

if ~isempty(col_seeds)&&(size(col_seeds,1)~=size(normed_data.Ec,2))
  error('Invalid column seed length');
end

if thr_row(1)<0||thr_col(1)<0
  warning('Negative thresholds. Are you sure about this?');
end

% Maybe add an error check for direction too

if strcmp(convergence,'cor')
  if cor_limit<=0
    warning('Non-positive correlation limit for convergence');
  end
end

if strcmp(convergence,'eps')  
  if eps>=1
    warning('`eps` limit for convergence greater than one');
  end
end

nSeeds = size(row_seeds,2) + size(col_seeds,2);

origTG = thr_row; %gene
origTC = thr_col;

if length(thr_row)~=1 && length(thr_row)~=nSeeds
  error('`thr_row` does not have the right length');
end
if length(thr_col)~=1 && length(thr_col)~=nSeeds
  error('`thr_col` does not have the right length');
end

thr_row = repmat(thr_row, nSeeds,1);
thr_col = repmat(thr_col, nSeeds,1);

% Putting the seeds together
%allSeeds = zeros(size(normed_data.Ec,1),1);
allSeeds = [];
if ~isempty(row_seeds)
  allSeeds = horzcat(allSeeds, row_seeds);
end
if ~isempty(col_seeds)
  col_seeds = isaRowfromCol(normed_data, col_seeds, thr_row(end-size(col_seeds,2)+1:end), direction{2}); %Last n= col size of col_seeds elements
  allSeeds = horzcat(allSeeds, col_seeds);
end

% All the data
rundata = struct('direction',{direction},'eps',eps,'cor_limit',...
                cor_limit,'maxiter',maxiter,'N',nSeeds,'convergence',...
                convergence,'prenormalize',normed_data.prenormalize,...
                'hasNaN',normed_data.hasNaN,'corx',corx,'unique',false);

% All the seed data; updated later
% seed data - matrix
%  --------------------------------------------------------------
%  | iterations | oscillations | thr_row | thr_col | freq | rob |
%  --------------------------------------------------------------
% seeddata = struct('iterations',NaN,'oscillation',0,'thr_row',thr_row, ...
%               'thr_col',thr_col,'freq',ones(nSeeds,1),'rob',NaN(nSeeds,1));

seeddata(:,3) = thr_row;
seeddata(:,4) = thr_col;
seeddata(:,1) = NaN;
seeddata(:,2) = 0;
seeddata(:,5) = ones(nSeeds,1);
seeddata(:,6) = NaN(nSeeds,1);

if isempty(allSeeds)
  isaresult = struct('Rows',allSeeds, 'Columns',...
                    zeros(size(normed_data.Ec,2),0), 'rundata', rundata,...
                    'seeddata',seeddata);
                  return;
end

%==========================================================================
% Choose convergence checking function %
if strcmp(convergence,'eps')
  checkConvergence = @epsConvergence; 
elseif strcmp(convergence,'cor')
  checkConvergence = @corConvergence;
elseif strcmp(convergence,'corx')
  if corx<2, error('Invalid `corx` value, should be atleast 2'); end
  rows_old = {};
  cols_old = {};
  idx_old = 1:corx;
  checkConvergence = @corxConvergence;
end

function res = epsConvergence(row_old, row_new, col_old, col_new)
tempr = row_old - row_new;
tempc = col_old - col_new;
for i=1:size(tempr,1)
  res(i) = all(abs(tempr(:,i))<eps) & all(abs(tempc(:,i))<eps);
end
res = res & ~isnan(res);
end

function res = corConvergence(row_old, row_new, col_old, col_new)
g_o = scale(row_old);
g_n = scale(row_new);
c_o = scale(col_old);
c_n = scale(col_new);

res = (sum(g_o*g_n,1)./(size(g_o,1)-1) > cor_limit) & ...
      (sum(c_o*c_n,1)./(size(c_o,1)-1) > cor_limit);
res = res & ~isnan(res);
end

function res = corxConvergence(row_old, row_new, col_old, col_new)
if iter<corx+1
  rows_old{end+1} = row_new;
  cols_old{end+1} = col_new;
  res = false(1,size(row_old,2));
else
  row_new = scale(row_new);
  col_new = scale(col_new);
  res = sum(rows_old{idx_old(1)}.*row_new,1)./(size(row_new,1)-1)>cor_limit...
      & sum(cols_old{idx_old(1)}.*col_new,1)./(size(col_new,1)-1)>cor_limit;
  rows_old{idx_old(1)} = row_new;
  cols_old{idx_old(1)} = col_new;
  idx_old = circshift(idx_old(:),-1)';
  res = res & ~isnan(res);
end
end
%==========================================================================

% Initializing
iter = 0;
index = 1:size(allSeeds,2);
row_old = allSeeds;
col_old = NaN(size(normed_data.Ec,2),nSeeds);
row_res = NaN(size(normed_data.Ec,1),nSeeds); % Change it zeros
col_res = NaN(size(normed_data.Ec,2),nSeeds); % Change it zeros

%--------------------------------------------------------------------------
% Main loop

while(true)
  
  iter = iter+1;
  [row_new col_new] = isaStep(normed_data, row_old, thr_row, thr_col, direction);

  % Mark converged seeds
  conv = checkConvergence(row_old, row_new, col_old, col_new); 
  
  % Mark all zero seeds
  zero0 = all(row_old); %tests whether all the elements are zero in column
  
%   % Mark oscillating ones, if requested
%   if oscillation && iter>1
%     new_fire = sum(row_new,1);
%   end
  osc = false;
  
  % these will be thrown out
  drop = conv|zero0|osc;
  
  % Drop the seeds to be dropped
  if(sum(drop)~=0)
    row_res(:,index(drop)) = row_new(:,drop);
    col_res(:,index(drop)) = col_new(:,drop);
%     seeddata.iterations(index(drop)) = iter;
    seeddata(index(drop),1) =iter; % iterations
    row_new = row_new(:,~drop); % Could be an error here
    col_new = col_new(:,~drop);
    thr_row = thr_row(~drop);
    thr_col = thr_col(~drop);
    if strcmp(convergence,'corx')
      rows_old = cellfun(@(x) x(:,~drop),rows_old, 'UniformOutput', false); % Another here
      cols_old = cellfun(@(x) x(:,~drop),cols_old, 'UniformOutput', false);
    end
    index = index(~drop);
  end
  
  if (size(row_new,2)==0||iter>maxiter), break; end
  
  row_old = row_new;
  col_old = col_new;
    
end % while ends
  
isaresult = struct('Rows',row_res, 'Columns',col_res, ...
                  'rundata',rundata, 'seeddata',seeddata);
end

function row_new = isaRowfromCol(normed_data, col_seeds, thr_row, direction)
Ec = normed_data.Ec;

if normed_data.hasNaN
  tempEc = Ec;
  tempEc(isnan(tempEc)) = 0;
  row_new = isafilter(tempEc*col_seeds, thr_row, direction);
else
  row_new = isafilter(Ec*col_seeds,     thr_row, direction);
end
end

