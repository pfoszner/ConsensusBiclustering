function parallelCoordinates(x, biClustResult, ClustNo, plotBoth, ...
                plotcol, compare, info, bothlab, order, ylab)
%Represent gene expression through gene and/or condition profile in a bicluster as lines

if nargin<10
  ylab = 'Value';
end

if nargin<9
  order = false;
end

if nargin<8
  bothlab = {'Rows', 'Columns'};
end

if nargin<7
  info = false;
end

if nargin<6
  compare = true;
end

if nargin<5
  plotcol = true;
end

if nargin<4
  plotBoth = false;
end

bicRows = biClustResult.Clust(ClustNo).rows;
bicCols = biClustResult.Clust(ClustNo).cols;

if order
  bicRows = bicRows(sortrows(sum(x(bicRows,bicCols),2)));
  bicCols = bicCols(sortrows(sum(x(bicRows,bicCols),1))); % Errror
end

if plotBoth
  figure;
  if compare
    mat = x(bicRows,:);
    subplot(2,1,1);
    plot(mat, 'Color','red');
    xlabel(bothlab{1});
    ylabel(ylab);
    hold all;
    plot(mat(:,bicCols));
    
    mat = x(:,bicCols)';
    subplot(2,1,2);
    plot(mat, 'Color','red');
    xlabel(bothlab{2});
    ylabel(ylab);
    hold all;
    plot(mat(:,bicRows));
  else
    subplot(2,1,1);
    plot(x(bicRows,bicCols));
    xlabel(bothlab{1});
    ylabel(ylab);
    
    subplot(1,1,2);
    plot(x(bicRows,bicCols)');
    xlabel(bothlab{1});
    ylabel(ylab);
  end
else
  if plotcol
    if compare
      mat = x(bicRows,:);
      plot(mat, 'Color','red');
      xlabel(bothlab{1});
      ylabel(ylab);
      hold all;
      plot(mat(:,bicCols));
    else
       plot(x(bicRows,bicCols));
       xlabel(bothlab{1});
       ylabel(ylab);
    end
  else
    if compare
      mat = x(:,bicCols)';
      subplot(2,1,2);
      plot(mat, 'Color','red');
      xlabel(bothlab{2});
      ylabel(ylab);
      hold all;
      plot(mat(:,bicRows));
    else
      plot(x(bicRows,bicCols)');
      xlabel(bothlab{1});
      ylabel(ylab);
    end
  end
end
if info
  title({'Bicluster', ['(rows=' num2str(length(bicRows)) '; ' 'columns=' num2str(length(bicCols)) ')']});
end
end