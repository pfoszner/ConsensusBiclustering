function biClustResult = plaid(matrix, cluster, fit_model, background,...
  background_layer, background_df, row_release, col_release, shuffle,...
   max_layers, iter_startup, iter_layer, verbose)
%Modified Plaid Algorithm by
% Heather Turner, Trevor Bailey, Wojtek Krzanowski,
% Improved biclustering of microarray data demonstrated through systematic performance tests,
% Computational Statistics & Data Analysis, Volume 48, Issue 2, 1 February 2005, Pages 235-254, ISSN 0167-9473, 10_1016/j_csda_2004_02_003_
% (http://www_sciencedirect_com/science/article/pii/S0167947304000295)
%
% Usage
% >> biClustResult =  bimaxBiclust(mat, alpha, delta, numBiclust, rand)
%
% Inputs
%   matrix                    - the data matrix
%   cluster                   - 'r', 'c' or 'b', to cluster rows, columns
%                                or both
%   fit_model                 - Model to fit each layer, usually linear.
%                               Default ['m' 'a' 'b'] => y ~ m+a+b
%                               Here  m: constant for all elements 
%                                     a: constant for all rows in bicluster
%                                     b: constant for all cols in bicluster
%   background                - if `true` the method will consider that the
%                               background layer (constant for all rows and
%                               columns) is present in the data matrix
%   background_layer          - if background=`true`, user specified
%                               background layer (matrix of dimension 
%                               `matrix`)
%   background_df             - degrees of freedom of background layer if
%                               background_layer is specified
%   shuffle                   - before a layer is added, its statistical
%                               significance is compared against the number
%                               of layers (obtained by random) defined by
%                               this parameter.
%   iter_startup              - number of iterations to find starting value
%   iter_layer                - number of iterations to find each layer
%   
%   row_release               - 0<scalar<1. Recommended:[0.5-0.7].
%                               Threshold to prune rows in the layers
%                               depending on row homogenity.
%   col_release               - same as row_release, with columns.
%   max_layers                - maximum number of layers to include in the
%                               model
%   verbose                   - if `true`, print extra information as the
%                               code progresses
%
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] if Row i is
%                     in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] if Col j is 
%                     in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices_
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India


if nargin<1
  error('input :  No matrix as input');
end

if nargin<2
  cluster          = 'b';
  fit_model        = ['m' 'a' 'b'];
  background       = true;
  background_layer = NaN;
  background_df    = 1;
  row_release      = 0.7;
  col_release      = 0.7;
  shuffle          = 3;
  
  max_layers       = 20;
  iter_startup     = 5;
  iter_layer       = 10;
  verbose          = true;
end

fix_layers          = [];
start_method        = 'convert';
iter_supervised     = [];
row_classes         = [];
col_classes         = [];
search_model        = [];
Z=matrix;

if (isempty(iter_startup)||isempty(iter_layer))
  error('Values for iter_startup and iter_layer have not been specified');
end

% Sorting out the input Z
if ismatrix(Z)
  Z(:,:,1) = Z;
end

[n,p,t] = size(Z);

if (isempty(row_classes)&&isempty(col_classes))
  iter_supervised = 0;
end

if (isempty(search_model))
  search_model = fit_model;
end

% Requirements on number of layers
if ~isempty(fix_layers)
  shuffle = 0;
  length = fix_layers+background;
else
  length = max_layers+background;
end

%##########################################################################
% Setting up variables to hold results
%##########################################################################

SS = NaN(length,1);
layer_df = NaN(length,1);
rows_released = NaN(length,1);
cols_released = NaN(length,1);
% A matrix with first column true and all other columns false
r = false(n,length);
r(:,1) = true;
% A matrix with first column true and all other columns false
k = false(p,length);
k(:,1) = true;
fits = cell(length,1);

%##########################################################################
% Fitting model begins
%##########################################################################
% background layer
if background
  if (~isnan(background_layer))
    fits{1} = background_layer;
    Z = bsxfun(@minus,Z,fits{1});
    temp = fits{1}.^2;
    SS(1) = sum(temp(:));
    if verbose
      disp(fprintf('layer: 0 \n',num2str(SS(1)),'\n'));
    end
    layer_df(1) = background_df;
  else
    fits{1} = fitLayer(Z, r(:,1),k(:,1), fit_model);
    Z = bsxfun(@minus, Z, fits{1});
    temp = fits{1}.^2;
    SS(1) = sum(temp(:));
    layer_df(1) = 1 + ismember('a',fit_model)*(n-1) + ismember('b',fit_model)*(p-1) + ismember('c',fit_model)*(t-1);
    
    if verbose
      disp(fprintf('layer: 0 \n',num2str(SS(1)),'\n'));
    end
  end
end

layer = +background;

%##########################################################################
% Bicluster Layers
%##########################################################################
count = 2;
while (layer < nanmin([fix_layers max_layers(:)])+background)
  if verbose
    disp(fprintf('layer: %d\n',layer));
  end
  [SS(count), r(:,count), k(:,count), fits{count}, layer_df(count), status(count), rows_released(count), cols_released(count)] = updatePlaid(Z,n,p,t,row_classes,col_classes,cluster,fit_model,search_model,row_release,col_release,shuffle,start_method,iter_startup,iter_layer,iter_supervised,verbose);
  
  % Stop if no clusters
  if SS(count)==0
    break;
  end
  count = count+1;
  % Else extract results; calculate new residual matrix
  layer = layer+1;
  Z(r(:,layer),k(:,layer),:) = bsxfun(@minus,Z(r(:,layer),k(:,layer),:),fits{layer});
end

%##########################################################################
% Results
%##########################################################################
if ~isempty(fix_layers)
  layer = fix_layers+background;
end

if verbose
  if layer==background
    disp('No biclusters found');
  else
    nr = sum(r(:,1:layer),1);
    nk = sum(k(:,1:layer),1);
    (fprintf('rows  cols df SS  MS Convergence  rows released  cols released')); 
    horzcat(nr',nk',layer_df(1:layer),SS(1:layer),SS(1:layer)./layer_df(1:layer), status(1:layer)', rows_released(1:layer), cols_released(1:layer))
  end
end

if layer<=1
  biClustResult = struct('RowxNum',NaN, 'NumxCol',NaN, ...
    'ClusterNo', 0, 'SS',0,'MS',0, 'Clust', NaN);
elseif layer==2
  RowxNum = r(:,2:layer);
  NumxCol = k(:,2:layer)';
  for i=1:1
    rows = find(RowxNum(:,i)>0);
    cols =find(NumxCol(i,:)>0);
    clustr(i) = struct('rows', rows, 'cols', cols);
  end
  biClustResult = struct('RowxNum',RowxNum,'NumxCol',NumxCol,'ClusterNo',1,'SS',SS(1:2),'MS',SS(1:2)./layer_df(1:2),'Clust',clustr);
else
  RowxNum = r(:,2:layer);
  NumxCol = k(:,2:layer)';
  for i=1:layer-1
    rows = find(RowxNum(:,i)>0);
    cols =find(NumxCol(i,:)>0);
    clustr(i) = struct('rows', rows, 'cols', cols);
  end
  biClustResult = struct('RowxNum',RowxNum,'NumxCol',NumxCol,'ClusterNo',layer-1,'SS',SS(1:layer),'MS',SS(1:layer)./layer_df(1:layer),'Clust',clustr);
end
end


%#########################################################################%
                        %% Helper functions %%
%#########################################################################%


function [SS, r, k, fits, layer_df, status, rows_released, cols_released]...
  = updatePlaid(Z, n, p, t, row_classes, col_classes, cluster, fit_model,...
                search_model, row_release, col_release, shuffle, ...
                start_method, iter_startup, iter_layer, iter_supervised,...
                verbose)
%% fit single bicluster

% set number of release iterations equal to number of layer iterations
if (~isempty(row_release)||~isempty(col_release))
  extra = round(iter_layer/2)*2;
else
  extra = 0;
end
flag = 1;
% set up the objects required
cluster_SS = zeros(shuffle+1,1);
status = 0;

for i=1:(shuffle+1)
  a = zeros(n,1);
  b = zeros(p,1);
  c = zeros(t,1);
  r_check = {1,1};
  k_check = {1,1};
  model = search_model;
  
  if i>1
    % permute genes and samples, within each time point separately 
   for ind=1:size(Z,3)
     g = Z(:);
     Z(:,:,ind) = reshape(g(randsample(n*p,n*p)),size(Z,1),size(Z,2));
   end
  end
  
  % This could be extended to add more methods like average etc.
  %# Rows
  if ismember(cluster, ['r' 'b'])
    if strcmp(start_method, 'convert')
      for ind=1:size(Z,3)
        r(:,ind) = mean(kmeansStart(Z(:,:,ind), iter_startup),2) >= 0.5; %row mean
      end
      if ~isempty(row_classes)
        temp_r = hist(r.*row_classes,max(row_classes)) >= row_grouped.*0.5;
        temp_r = temp_r(row_classes);
        if sum(temp_r(:))~=0
          r = temp_r;
        else 
          row_classes = [];
          if i==1
            fprintf('Row starting values all converted to zero.\nReverting to unsupervised iterations.');
          end
        end
      end
      if sum(r(:))>n/2
        r = ~r;
      end
    end
  else
    r = true(n,1);
  end
  
  % This could be extended to add more methods like average etc.
  %# Columns
  if ismember(cluster, ['c' 'b'])
    if strcmp(start_method, 'convert')
      for ind=1:size(Z,3)
        k = mean(kmeansStart(Z(:,:,ind), iter_startup,true),2) >= 0.5; %row mean
      end
      if ~isempty(col_classes)
        temp_k = hist(k.*col_classes,max(col_classes)) >= col_grouped.*0.5;
        temp_k = temp_k(col_classes);
        if sum(temp_k(:)~=0)
          k = temp_k;
        else
          col_classes = [];
          if i==1
            disp('Column starting values all convertes to zero. Reverting to unsupervised iterations.');
          end
        end
      end
      if sum(k(:))>p/2
        k = ~k;
      end
    end
  else
    k = true(p,1);
  end

  %------------------------------------------------------------------------
  j = 0;
  while(j<=iter_layer+extra)
    if (0<j && j<=iter_layer)
      % update cluster membership parameters, r 
      if ismember(cluster, ['r' 'b'])
        if (j<=iter_supervised && ~isempty(row_classes))
          % supervised updates
          r = sum(sum((Z(:,k,:)-makeLayer(m,a,b(k),c)).^2,1),2) < ...
              sum(sum(Z(:,k,:).^2,1),2); % Rowsum
          r = r(row_classes);
          if sum(r(:) == 0)
            row_classes = [];
            fprintf('Layer %d: no rows clustered - reverting to unsupervised', i);
          end
        end
        if (j>iter_supervised || isempty(row_classes))
          r = sum((Z(:,k,:)-makeLayer(m,a,b(k),c)).^2,2) < sum(Z(:,k,:).^2,2); %Rowsum
        end
      end
      n2 = sum(r(:));
      if n2==0
        break;
      end
      % update cluster membership parameters, k
      if ismember(cluster, ['c' 'b'])
        if (j<=iter_supervised && ~isempty(col_classes))
          % supervised updates
          k = sum(sum(permute((Z(r,:,:)-makeLayer(m,a(r),b,c)).^2, [2 1 3]),2),2) < sum(sum(permute(Z(r,:,:).^2,[2 1 3]),2),2); %Rowsum
          k = k(col_classes);
          if sum(k(:)==0)
            col_classes=[];
            fprintf('Layer %d: no columns clustered - reverting to unsupervised',i);
          end
        end
        if (j>iter_supervised || isempty(col_classes))
          k = sum(permute((Z(r,:,:)-makeLayer(m,a(r),b,c)).^2,[2 1 3]),2) < sum(permute(Z(r,:,:).^2,[2 1 3]),2);
        end
      end
    end
    
    if (j>=iter_layer+1 && ~isempty(row_release) && mod(j-iter_layer,2)==1)
      % row release
      if resdf==0
        r = false(n,1);
      else
        r(r) = (1/resdf) * sum((Z(r,k,:)-makeLayer(m,a(r),b(k),c)).^2,2)...
          < (1-row_release)/totdf * sum(Z(r,k,:).^2,2); % 
      end
    end
    n2 = sum(r(:));

    if (j>=iter_layer+1 && ~isempty(col_release) && mod(j-iter_layer,2)==0)
      % column release
      if (totdf==0||resdf==0)
        k = false(p,1);
      else
        k(k) = (1/resdf) * sum(sum((Z(r,k,:)-makeLayer(m,a(r),b(k),c)).^2,1),2) < (1-col_release)/totdf * sum(sum(Z(r,k,:).^2),2); %Recheck for bug
      end
    end
    p2 = sum(k(:));
    
    if (n2==0 || p2==0)
      if i==1
        if verbose
          disp([j sum(r(:)), sum(k(:))]);
        end
        stopnow = true;
        n_iter = j+1;
      end
      break;
    end
    
    
    % Skipping to last iteration if already converged
    if (j >= iter_supervised && j <= iter_layer)||(j > iter_layer && mod(j-iter_layer,2)==0)
      tempr = 1:n;
      tempk = 1:p;
      if flag==1
        r_check{2} = tempr(r); 
        k_check{2} = tempk(k); 
        flag=2;
      else
        r_check{1} = tempr(r);
        k_check{1} = tempk(k);
        flag=1;
      end
      if (isequal(size(r_check{1},2), size(r_check{2},2)))&&isequal(size(k_check{1},2),size(k_check{2},2))
        if (isequal(r_check{1}, r_check{2}))&&isequal(k_check{1}, k_check{2})
          if i==1 && j<=iter_layer
            n_iter = j;
            status = i;
          end
          if (j>iter_layer)
            j = iter_layer+extra;
          else
            j = iter_layer;
          end
        end
      end
    end
    
    if j==iter_layer
      % use fit_model for final model/basis of row & col
      model = fit_model;
      % save no. of rows & cols in order to calculate no. released
      n_rows = n2;
      n_cols = p2;
    end
    % calculate d.f. for row/col release
    if j>=iter_layer
      totdf = n2*p2*t;
      resdf = totdf - (1+ismember('a',model)*(n2-1)+ismember('b', model)*(p2-1)+ismember('c',model)*(t-1));
      if i==1 && resdf==0 && verbose
        disp('Zero residual degrees of freedom.');
      end
    end
    
    % update layer effects
    m = mean2(Z(r,k,:));
    if any(ismember('a', model))
      a(r) = mean(bsxfun(@minus, Z(r,k,:),m),2); %Row mean 
      a(~r) = 0;
    end
    if any(ismember('b', model))
      b(k) = mean(bsxfun(@minus,Z(r,k,:),m),1)'; %Column mean 
      b(~k) = 0;
    end
    if any(ismember('c', model))
      c = mean(bsxfun(@minus,Z(r,k,:),m),1); %Column mean
    end
    
    if i==1 && verbose
      disp([j sum(r(:)), sum(k(:))]);
    end
    j = j+1;
  end
  
  if n2==0||p2==0
    cluster_SS(i) = 0; % Maybe NaN
  else
    cluster_SS(i) = sumsqr(makeLayer(m,a(r),b(k),c));
  end
  if i==1
    % Save results for candidate layers
    if ~exist('n_iter'),    n_iter = j-1;    end
    if verbose,      disp(n_iter);    end
    
    id = {squeeze(r), squeeze(k)};
    if n2==0||p2==0
      fits = NaN;
      layer_df = NaN;
      rows_released = NaN;
      cols_released = NaN;
    else
      fits = makeLayer(m, a(r),b(k),c);
      layer_df = 1 + ismember('a', model)*(n2-1) + ismember('b', model)*(p2-1) + ismember('c', model)*(t-1);
      if ~isempty(row_release)||cluster=='c'
        rows_released = n_rows - n2;
      else
        rows_released = NaN;
      end
      if ~isempty(col_release)||cluster=='r'
        cols_released = n_cols - p2;
      else
        cols_released = NaN;
      end
    end
  end
  if exist('stopnow')
    break;
  end
end

if verbose
  disp(cluster_SS);
end

if shuffle>0
  if cluster_SS(1)>max(cluster_SS(2:end))
    cluster_SS = cluster_SS(1);
  else
    cluster_SS = 0;
  end
end
SS = cluster_SS;
r = id{1};
k = id{2};

end

%##########################################################################

function x = kmeansStart(Z, iter_startup, transpose)
% find k-means starting values for rows (cols if transpose = true)
if nargin<3
  transpose = false;
end
try
  if transpose
    km = kmeans(Z',2,'options',statset('MaxIter',iter_startup));
  else
    km = kmeans(Z,2,'options',statset('MaxIter',iter_startup));
  end
  km_size = histc(km,unique(km));
  for i=1:length(km)
    if km_size(1)<km_size(2)
      if km(i)==1
        x(i,1) = true;
      else
        x(i,1) =false;
      end
    else
      if km(i)==2
        x(i,1) = true;
      else
        x(i,1) = false;
      end
    end
  end
catch err
  x = false(size(Z,1+transpose),1);
end
end

%##########################################################################

function out = fitLayer(Z, r, k, model)
% Function to fit layer, given residuals from all other layers in the model

Z = Z(r,k,:);
m = mean2(Z); % Mean of entire matrix

if (any(ismember('a', model)))
  a = mean(bsxfun(@minus,Z,m),2); % row means
else
  a = zeros(size(Z,1),1);
end

if (any(ismember('b', model)))
  b = mean(bsxfun(@minus,Z,m),1)'; % Column Means
else
  b = zeros(size(Z,2),1);
end

if (any(ismember('c', model)))
  c = mean(bsxfun(@minus,Z,m),1); % Column means
else
  c = zeros(size(Z,3),1);
end
out = makeLayer(m,a,b,c);
end

%##########################################################################

function out = makeLayer(m,a,b,c)
% Function to construct layer given the fitted effects

ma = bsxfun(@plus,m,a);
out = zeros(length(ma),length(b),length(c));
for i=1:size(ma)
  for j=1:size(b)
    for k=1:size(c)
      out(i,j,k) = ma(i)+b(j)+c(k);
    end
  end
end
end
