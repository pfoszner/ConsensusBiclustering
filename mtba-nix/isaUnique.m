function isaresult = isaUnique(normed_data, isaresult, method,...
                               ignore_div, cor_limit, neg_cor, drop_zero)
% Create a unique set of biclusters, by removing all biclusters that are
% similar to others.
%
% Usage: isaresult = isaUnique(normed_data, isaresult, method,...
%                              ignore_div, cor_limit, neg_cor, drop_zero)
%
% Inputs
%   normed_data   -   the normalized data consisting of two matrices, first
%                     one transposed. Also contains information if the data
%                     was prenormalized or contains NaN.
%   isaresult     -   the result of an ISA run, which is a set of
%                     biclusters
%   method        -   String, giving the method to be used for getting
%                     similarities. Right now only 'cor' is implemented
%   ignore_div    -   if `true` then divergent biclusters will be removed.
%                     Default:true
%   cor_limit     -   Number giving the correlation limit for 'cor' method.
%                     Default:0.9
%   neg_cor       -   if `true`, then 'cor' method considers the absolute
%                     value of correlation
%   drop_zero     -   if `true` drops all biclusters which have all zero
%                     scores
% Output
%   isaresult: A structure consisting of
%     Rows          - the row components of biclusters, every column
%                     corresponding to a bicluster.
%     Columns       - the column components of biclusters, in the same
%                     format as rows.
%     rundata       - structure with information about ISA runs.
%     seeddata      - matrix with information about the biclusters.
%
% Because of the nature of ISA algorithm, isaIterate does not give a set of
% unique biclusters. So this method filters and removes the ones which are
% very similar.
%
% See Also: itersa
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin<3
  method = 'cor';
  ignore_div = true;
  cor_limit = 0.9;
  neg_cor = true;
  drop_zero = true;
end

if size(isaresult.Rows,2)==0, return; end

% Drop divergent seeds
if ignore_div
  invalid = isnan(isaresult.seeddata(:,1));
  if (any(invalid(:)))
    valid = ~invalid;
    isaresult.Rows = isaresult.Rows(:,valid);
    isaresult.Columns = isaresult.Columns(:,valid);
    isaresult.seeddata = isaresult.seeddata(valid,:);
  end
end

if size(isaresult.Rows,2)==0, return; end

% Drop zero seeds
if drop_zero
  valid = any(isaresult.Rows);
  if ~all(valid(:))
    isaresult.Rows = isaresult.Rows(:,valid);
    isaresult.Columns = isaresult.Columns(:,valid);
    isaresult.seeddata = isaresult.seeddata(valid,:);
  end
end
if size(isaresult.Rows,2)==0, return; end

if strcmp(method, 'cor')
  % We reorder the results to keep mpdules with higher (row,column)
  % threshold
  [~,ord] = sortrows(isaresult.seeddata(:,[3 4])); % order by first column, decide ties by second
%   ord = flipdim(ord,1); % Decreasing order
  isaresult.Rows = isaresult.Rows(:,ord);
  isaresult.Columns = isaresult.Columns(:,ord);
  isaresult.seeddata = isaresult.seeddata(ord,:);
  if neg_cor, ABS=@abs; else ABS=@(x) x; end
  %         cm(1,:) = min(ABS(corrcoef(isaresult.Rows))) ; % Error
  %         cm(2,:) = min(ABS(corrcoef(isaresult.Columns))); % Error
  cm = min(ABS(corrcoef(isaresult.Rows)), ABS(corrcoef(isaresult.Columns)));
  cm(find(tril(cm))) = 0;
  uni = all(bsxfun(@lt, cm, cor_limit),1);
  freq = zeros(size(cm,1),1); %Preallocating
  for i=1:size(cm,1)
    tempf = isaresult.seeddata(bsxfun(@ge,cm(i,:),cor_limit),5);
    freq(i,1) = sum(tempf(:)) + isaresult.seeddata(i,5);
  end
  %         freq = arrayfun(@(x) sum(sum(isaresult.seeddata((cm(x,:) >= cor_limit),5))), 1:size(cm,1))+isaresult.seeddata(:,5);
elseif strcmp(method,'round') %TODO
end

isaresult.Rows = isaresult.Rows(:,uni);
isaresult.Columns = isaresult.Columns(:,uni);
isaresult.seeddata = isaresult.seeddata(uni,:);
isaresult.seeddata(:,5) = freq(uni);
isaresult.rundata.unique = true;


end
