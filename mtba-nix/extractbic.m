function result = extractbic(resinfo,matrix,nRow,nCol,p)
%% Function to extract Biclusters from the output of infocontent
% 
% Input :
%       resinfo  -  result of infocontent(Structure Consisting of loading 
%                   matrix and factor matrix)
%        matrix  -  matrix obtained after preprocessing
%          nRow  -  number of rows in the matrix
%          nCol  -  number of coloums in the matrix
%             p  -  number of biclusters to be extracted
%             
% Output :
%        result  -  final bicluster model 
% 
% Author: Shruti jain, 2014
%        
% Contact: sjain@iitk.ac.in, sweetushruti963@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India
%
%% Extracting biclusters
noL = resinfo.noL;
nZ = resinfo.nZ;
thresZ = 0.5;
thresL = 0;

for i = 1:p
    thresL = thresL+(sum(noL(:,i).^2))*(sum(nZ(i,:).^2));
end
thresL = thresL/(nRow*nCol*p);
thresL = sqrt(thresL)/thresZ;

for i = 1:p
    model(i) = struct('GeneIndex',[],'TissueIndexp',[],'TissueIndexn',[]);
end

for i = 1:p
    gene = find(abs(noL(:,i))>thresL);
    model(i).GeneIndex = gene;
    samplep = find(nZ(i,:)>thresZ);
    samplen = find(nZ(i,:)<-thresZ);
    sump = sum(abs(nZ(i,samplep)));
    sumn = sum(abs(nZ(i,samplen)));
    if sump>=sumn
        model(i).TissueIndexp = samplep;
        model(i).TissueIndexn = samplen;
    else
        model(i).TissueIndexp = samplen;
        model(i).TissueIndexn = samplep;
    end
end

result = model;

end        