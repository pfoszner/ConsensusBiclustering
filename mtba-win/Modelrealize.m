function finalmodel = Modelrealize(model,dataset,ara,t,InT,InG)

if size(model,2) == 0
    finalmodel = [];
end

best = 0;
pBest = -1;

for i=2:size(model,2)
    if model(i).p > pBest
        best = i;
        pBest = model(i).p;
        isInT = InT(i-1,:);
        isInG = InG(i-1,:);
        arrayTissue = ara{1,i-1};
    end
end

finalmodel = extend(dataset,t,model(best),arrayTissue,isInT,isInG);
end