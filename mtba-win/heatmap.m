function Sortedmatrix = heatmap(matrix, biCluster, ClusterNo)
% Method to display a heatmap after sorting the input
%
% Usage
% >> Sortedmatrix = heatmap(matrix, biCluster, ClusterNo)
%
% Inputs:
%   matrix                  - The data matrix where the bicluster is to
%                             be drawn.
%   biCluster (optional)    - A bicluster result set. If absent the data
%                             matrix is drawn as a heatmap, without any
%                             reordering.
%   ClusterNo               - Index of cluster to be plotted
%
% Outputs:
%   Sortedmatrix    - Sorted matrix
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

if nargin < 1
    error('input :  No matrix as input');
end

if nargin < 2
    Sortedmatrix = matrix;
    figure;
    colormap(winter);
    imagesc(matrix);
    set(gca,'FontSize',14,'FontWeight','bold');
end

if nargin == 3
    if (biCluster.ClusterNo<ClusterNo || ClusterNo<0)
        error('Wrong Cluster Number: Cluster does not exist');
    end
    
    % Sorting Matrix
    % According to the Columns
    [~,soc] = sort(biCluster.NumxCol(ClusterNo,:),'descend');
    socmatrix = matrix(:,soc);
    % According to the Rows
    [~,sor] = sort(biCluster.RowxNum(:,ClusterNo),'descend');
    Sortedmatrix = socmatrix(sor,:);
    % Plotting
    figure;
    colormap(winter);
    imagesc(Sortedmatrix);
    set(gca,'FontSize',14,'FontWeight','bold');
    hold on;
    plot([length(biCluster.Clust(ClusterNo).cols)+0.5, ...
        length(biCluster.Clust(ClusterNo).cols) + 0.5],...
        [0.5, length(biCluster.Clust(ClusterNo).rows)+0.5], 'Color', 'r');
    plot([0.5, length(biCluster.Clust(ClusterNo).cols)+0.5],...
        [length(biCluster.Clust(ClusterNo).rows)+0.5, ...
        length(biCluster.Clust(ClusterNo).rows)+0.5],'Color','r');
    hold off;
end

if nargin > 3
    error('input: Wrong number of input');
end

