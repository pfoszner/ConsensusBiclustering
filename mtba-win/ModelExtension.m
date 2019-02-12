function extendedModel = ModelExtension(parent,t,dataset,arrayTissue,isInG)

a=parent.a;
b=parent.b;
p = getScoreForExtension(dataset,t,arrayTissue,a,b,isInG);
extendedModel = struct('p',p);
end

function p = getScoreForExtension(dataset,t,arrayTissue,a,b,isInG)

noRow = size(dataset,1);
noCol = size(dataset,2);
[so,ind] = sort(dataset,2);
[s1,D] = sort(ind,2);
ta = arrayTissue(a);
tb = arrayTissue(length(arrayTissue)-b+1);

if a+b==length(arrayTissue)-1
    k=0;
    for row=1:noRow
        if (isInG(row) && D(row,ta)<D(row,t) && D(row,t)<D(row,tb))
            k=k+1;
        end
    end
    p=k/noRow;
else
    commonDenominator = nchoosek(noCol,length(arrayTissue));
    B = 1/(nchoosek(noCol,a+b+1)*factorial(a+b+1));
    for i=1:noRow
        A(i) = 0;
        if (isInG(i) && D(i,ta)<D(i,t) && D(i,t)<D(i,tb))
            gi = 0;
            if a==b
                gi = D(i,tb) - D(i,t) - 1;
            else
                gi = D(i,t) - D(i,ta) - 1;
            end
            if gi >= length(arrayTissue)-a-b-1
                A(i) = nchoosek(gi,length(arrayTissue)-a-b-1)/commonDenominator;
            else
                A(i)=0;
            end
        end
            diff(i) = A(i) - B;
    end
        
        p = 0.05;              % Starting Point for Newton Iteration
        for iteration = 1:20
            f = -noRow;
            dfdp = 0;
            for i = 1:noRow
                denominator(i) = diff(i)*p + B;
                if denominator(i)==0 && A(i)~=0
                    disp('WARNING: Division by 0 in Newton Iteration.');
                    f = 0;
                    dfdp = 1;
                else if A(i)~=0
                        f = f + (A(i)/denominator(i));
                        dfdp = dfdp - (A(i)*diff(i)/(denominator(i)^2));
                    end
                end
            end

            if dfdp==0
                p = 0;
            else
                q = p;
                p = p - (f/dfdp);
                if abs(f/dfdp) < 0.0001
                    break;
                end
                if p<0 && iteration<20
                    p = q/2;
                end
            end
        end

        if p<0 || p>1
            p=0;
        end

    end
                
    
end

function extend(t,a,b,arrayTissue,isInT)

if a+b==length(arrayTissue)
    error('Complete Models cannot be extended!');
end

if a==b
    arrayTissue(a+1) = t; a=a+1;
else
    arrayTissue(length(arrayTissue)-b) = t; b = b+1;
end

isInT(t) = 1;
updateCompatibilityInfo()
end