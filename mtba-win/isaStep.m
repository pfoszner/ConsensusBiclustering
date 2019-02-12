function [row_new col_new]  = isaStep(normed_data, rows, thr_row, thr_col, direction)
% Step filter the normalized data
% (Internal function for ISA algorithm)
%
% See Also: isaIterate, isafilter
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

Ec = normed_data.Ec;
Er = normed_data.Er;

% Direction validity could be added


if normed_data.hasNaN
  tempEr = Er;
  tempEr(isnan(tempEr)) = 0;
  tempEc = Ec;
  tempEc(isnan(tempEc)) = 0;
  col_new = isafilter(tempEr*rows,    thr_col, direction{1});
  row_new = isafilter(tempEc*col_new, thr_row, direction{2});
else
  col_new = isafilter(Er*rows,    thr_col, direction{1});
  row_new = isafilter(Ec*col_new, thr_row, direction{2});
end
  
end