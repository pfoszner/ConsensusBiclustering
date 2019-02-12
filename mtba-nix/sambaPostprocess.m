function filteredBicluster = Postprocess(Bicluster,L)


for i=1:size(Bicluster,1)
    for j=(1:size(Bicluster{i,1},2))
        if i==1 && j==1
            combBic = struct('rows',Bicluster{i,1}(1,j).rows,'column',Bicluster{i,1}(1,j).column);
        else
            combBic = [combBic,Bicluster{i,1}(1,j)];
        end
    end  
end

j=2;
while j <= size(combBic,2)
    for j_1 = 1:size(combBic,2)
        if j_1~=j
            sharedCond = intersect(combBic(1,j).column,combBic(1,j_1).column);
            sharedGene = intersect(combBic(1,j).rows,combBic(1,j_1).rows);
            if (size(sharedGene,2)*size(sharedCond,2))/(size(combBic(1,j).column,2)*size(combBic(1,j).rows,2)) > (L/100)
                combBic = combBic([1:j-1, j+1:size(combBic,2)]);
                j=j-1;
                break;
            end
        end        
    end
    j=j+1;
end

filteredBicluster = combBic;
end