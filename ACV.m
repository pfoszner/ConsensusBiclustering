function value = ACV(Amatrix)

    [dim1, dim2] = size(Amatrix);
    
    RowValue = 0;

    div = 0;

    for i = 1:dim1 - 1
    
        iRow = Amatrix(i,:);
        iAverage = mean(iRow);
        denominator1 = sum((iRow - iAverage) .* (iRow - iAverage));

        denominator1 = sqrt(denominator1);

        for j = i + 1:dim1
        
            jRow = Amatrix(j,:);
            jAverage = mean(jRow);
            nominator = sum((iRow - iAverage) .* (jRow - jAverage));

            denominator2 = sum((jRow - jAverage) .* (jRow - jAverage));
            denominator = denominator1 * sqrt(denominator2);

            if (denominator ~= 0)
                RowValue = RowValue + abs(nominator / denominator);
            else
                RowValue = RowValue + 1;
            end

            div = div + 1;
        end
    end

    RowValue = RowValue / div;

    ColumnValue = 0;

    div = 0;

    for i = 1:dim2 - 1
    
        iCol = Amatrix(:,i);
        iAverage = mean(iCol);
        denominator1 = sum((iCol - iAverage) .* (iCol - iAverage));
        denominator1 = sqrt(denominator1);

        for j = i + 1:dim2
        
            jCol = Amatrix(:,j);

            jAverage = mean(jCol);

            nominator = sum((iCol - iAverage) .* (jCol - jAverage));

            denominator2 = sum((jCol - jAverage) .* (jCol - jAverage));
            denominator = denominator1 * sqrt(denominator2);

            if (denominator ~= 0)
                ColumnValue = ColumnValue + abs(nominator / denominator);
            else
                ColumnValue = ColumnValue + 1;
            end
            div = div + 1;
        end
    end

    ColumnValue = ColumnValue / div;
    
    if (RowValue > ColumnValue)
        value = RowValue;
    else
        value = ColumnValue;
    end

end