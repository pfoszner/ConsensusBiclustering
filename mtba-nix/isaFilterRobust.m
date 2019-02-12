function isares = isaFilterRobust(data, normed_data, isares, perms, ...
                                  row_seeds, col_seeds)
% Wrapper to robustness for ISA
%
% Usage
% >> isares = isaFilterRobust(data, normed_data, isares, perms, ...
%                                   row_seeds, col_seeds)
% Inputs
%   data            -   the original, non normalized data matrix
%   normed_data     -   the normalized input data, calculated with
%                       isaNormalize
%   isares          -   the result of ISA, coming from isaIterate
%   perms           -   number of permutations to perform on input data
%   row_seeds       -   Row seeds for ISA run on scambled data
%   col_seeds       -   Column seeds for ISA run on scrambled data
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
% See Also: itersa, generateSeeds, isaIterate, isarobustness, isaNormalize
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if perms<=0
  error('input: Number of permutations must be positive');
end

if isempty(isares.Rows)
  return;
end

if length(unique(isares.seeddata(:,3)))~=1 %thr_row
  warning('Different row thresholds, using only the first one');
end

if length(unique(isares.seeddata(:,4)))~=1 %thr_col
  warning('Different column thresholds, using only the first one');
end

isares.seeddata(:,6) = isarobustness(normed_data, isares.Rows, isares.Columns); %rob

if isempty(row_seeds)
  row_seeds = generateSeeds(size(isares.Rows,1),isares.rundata.N);
end

if isempty(col_seeds)
  col_seeds = zeros(size(isares.Columns,1),0);
end

rob_max = 0;

for i=1:perms
  eggs = data(randperm(numel(data)));
  data_scrambled = reshape(eggs, size(data));

  normed_data_scrambled = isaNormalize(data_scrambled);
  
  permres = isaIterate(normed_data_scrambled, row_seeds, col_seeds,...
                       isares.seeddata(i,3), ... % thr_row
                       isares.seeddata(i,4), ... % thr_col
                       isares.rundata.direction, ...
                       isares.rundata.maxiter, ...
                       isares.rundata.convergence, ...
                       isares.rundata.cor_limit, isares.rundata.eps, ...
                       isares.rundata.corx ...
                       ); 
  valid = any(permres.Rows~=0,1);
  valid = valid & any(permres.Columns~=0,1);
  
  permres.Rows = permres.Rows(:,valid);
  permres.Columns = permres.Columns(:,valid);
  permres.seeddata = permres.seeddata(valid,:);
  
  rob2 = isarobustness(normed_data_scrambled, permres.Rows, permres.Columns);
  rob_max = nanmax([rob2 rob_max]);
end

keep = isares.seeddata(:,6) > rob_max; %rob

isares.Rows = isares.Rows(:,keep);
isares.Columns = isares.Columns(:,keep);
isares.seeddata = isares.seeddata(keep,:);

if size(isares.seeddata,1)>0, isares.seeddata(:,7)= rob_max; end
isares.rundata.rob_perms = perms;
end  