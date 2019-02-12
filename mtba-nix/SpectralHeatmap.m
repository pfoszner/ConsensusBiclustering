function SpectralHeatmap(matrix,row_cluster,column_cluster)
% Function to plot Heat Map of Spectral Co-clustered data
% Inputs:
%        matrix           :       Input matrix [n*m] of co-occurences of n instances (rows) and m
%                                 features (columns).
%        row_cluster      :       Input vector containing cluster number of instances(rows).
%        column_cluster   :       Input vector containing cluster number of features(columns).
%
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

% Input Argument Check
if nargin < 1
    error('Input matrix; row and column cluster number vector not specified.');
end
if (isempty(matrix))
    error('Input matrix is empty.')
end
if (length(find(sum(abs(matrix),1)==0))>0)
    error('Input matrix has empty instances (zero values along one or more columns). Replace atleast one of them with non-zero values.');
end
if (length(find(sum(abs(matrix),1)==0))>0)
    error('Input matrix has empty features (zero values along one or more rows). Replace atleast one of them with non-zero values.')
end
if nargin < 2
    error('Row and column cluster numbers not specified.');
end
if nargin < 3
    error('Column cluster numbers not specified.');
end

%Sorting the column cluster according to their cluster number
sorted_column=sort(column_cluster);
x_index=[];
uni_c=unique(sorted_column);
for i=1:length(unique(sorted_column))
    x_index=[x_index ; find(column_cluster==uni_c(i))];
end

%Sorting the row cluster according to their cluster number
sorted_row=sort(row_cluster);
y_index=[];
uni_r=unique(sorted_row);
for i=1:length(unique(sorted_row))
    y_index=[y_index ; find(row_cluster==uni_r(i))];
end

%Plotting the heatmap
sorted_matrix=matrix(y_index,x_index);
figure; 
colormap(winter);
imagesc(sorted_matrix);
set(gca,'fontsize',14,'fontweight','bold');
hold on;

%Following code draws a rectangular boundary around the biclusters found
num=0;
for i=1:length(uni_c)
    if sum(uni_r==uni_c(i))
        num=num+1;
        %text(find(sorted_column==uni_c(i),1)-0.5  , find(sorted_row==uni_c(i),1)-0.5 + (sum(sorted_row==uni_c(i))/2),sprintf('Bi-cluster %d',num),'color','red','fontsize',18);
        rectangle('position',[find(sorted_column==uni_c(i),1)-0.5 find(sorted_row==uni_c(i),1)-0.5 sum(sorted_column==uni_c(i)) sum(sorted_row==uni_c(i))],'LineWidth',2,'EdgeColor','red');
    end
end
end
