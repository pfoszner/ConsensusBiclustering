function result = fabic(matrix,psi,L,lapla,cyc,alpha,spl,spz,scale,lap,nL,lL,bL,non_negative)
%% Function for updation of sparse loadings
%
% Inputs :
%             matrix  -  Input datamatrix
%                psi  -  noise variance
%              lapla  -  variational parameter
%                cyc  -  Number of Iterations;default = 500
%              alpha  -  sparseness loadings;default = 0.01
%                spl  -  sparseness prior loadings;default = 0;
%                spz  -  sparseness factor;default = 0.5
%              scale  -  loading vectors are scaled in each iteration to the given variance. 0.0 indicates non scaling; default = 0        
%                lap  -  minimal value of the variational parameter, default = 1
%                 nL  -  maximal number of biclusters at which a row element can participate; default = 0(no limit)
%                 lL  -  maximal number of row elements per bicluster; default = 0 (no limit)
%                 bL  -  cycle at which the nL or lL maximum starts; default = 0 (start at the beginning)
%       non_negative  -  non_negative loading and factors :0 for no, 1 for yes:default = 0
%  
% Output :
%             result  -  final bicluster model 
%
% See Also : infocontent.m, extractbic.m
%
% Author: Shruti jain, 2014
%        
% Contact: sjain@iitk.ac.in, sweetushruti963@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India
%


nRow = size(matrix,1);
nCol = size(matrix,2);
p = size(L,2);

%xx = zeros(1,nRow);
e_sx_n = zeros(1,p);
e_ssxx_n = zeros(1,p);
ichol = zeros(p,p);
sum2 = zeros(p,p);
%tLpsiL = zeros(p,p);
LpsiL = zeros(p,p);
Lpsi = zeros(nRow,p);
sum1 = zeros(nRow,p);

spl = -spl;
spz = -spz;
in = 1/nCol;

eps = 0.001;
epsl = 10^-10;

if lap<eps
    lap = eps;
end

xx = in*(sum(matrix.*matrix,2));
xx(logical(xx<eps)) = eps;

for cycle = 1:cyc
    
   Lpsi = bsxfun(@rdivide,L',psi')';
   for i1 = 1:p
       for i2 = i1:p
           LpsiL(i1,i2) = sum(L(:,i1).*L(:,i2)./psi);
           LpsiL(i2,i1) = sum(L(:,i1).*L(:,i2)./psi);
       end
   end
   
   sum1(:,:) = 0;
   sum2(:,:) = 0;
   sum2(logical(eye(p,p))) = eps;
   
   for j = 1:nCol
       tLpsiL = LpsiL;
       for i1 = 1:p
           tLpsiL(i1,i1) = LpsiL(i1,i1)+lapla(j,i1);
       end
       
       for i1 = 1:p
           for i2 = i1:p
              s = tLpsiL(i1,i2);
              for i3 = i1-1:-1:1
                  s = s-tLpsiL(i1,i3)*tLpsiL(i2,i3);
              end
              ichol(i1,i2) = 0;
              ichol(i2,i1) = 0;
              if i1==i2
                  tLpsiL(i1,i1) = sqrt(s);
              else
                  tLpsiL(i2,i1) = s/tLpsiL(i1,i1);
              end
           end
       end
       
       for i1 = 1:p
           for i2 = 1:i1
               if i1 == i2
                   s = 1;
               else
                   s = 0;
               end
               for i3 =i1-1:-1:i2
                   s = s-tLpsiL(i1,i3)*ichol(i2,i3);
               end
               ichol(i2,i1) = s/tLpsiL(i1,i1);
           end
       end
      
       for i1 = p:-1:1
           for i2 = 1:i1
               if i1<i2
                   s = 0;
               else
                   s = ichol(i2,i1);
               end
               for i3 = i1+1:p
                   s = s-tLpsiL(i3,i1)*ichol(i2,i3);
               end
               ichol(i1,i2) = s/tLpsiL(i1,i1);
               ichol(i2,i1) = s/tLpsiL(i1,i1);
           end
       end
       
       for i3 = 1:p
           e_ssxx_n(i3) = sum(Lpsi(:,i3).*matrix(:,j));
       end
       
       for i1 = 1:p
           t = sum(ichol(i1,:).*e_ssxx_n);
           if t<epsl && non_negative>0
               e_sx_n(i1) = 0;
           else
               e_sx_n(i1) = t;
               sum1(:,i1) = sum1(:,i1)+(t*matrix(:,j));
           end
       end
       
       for i1 = 1:p
           for i2 = 1:p
               s = ichol(i1,i2)+e_sx_n(i1)*e_sx_n(i2);
               sum2(i1,i2) = sum2(i1,i2)+s;
               if i1 == i2
                   e_ssxx_n(i1) = s;
               end
           end
       end
       
       for i1 = 1:p
           s = (eps+e_ssxx_n(i1))^spz;
           if s<lap
               lapla(j,i1) = lap;
           else
               lapla(j,i1) = s;
           end
       end
   end
   
   for i1 = 1:p
       for i2 = i1:p
           s = sum2(i1,i2);
           for i3 = i1-1:-1:1
               s = s-sum2(i1,i3)*sum2(i2,i3);
           end
           ichol(i1,i2) = 0;
           ichol(i2,i1) = 0;
           if i1 == i2
               sum2(i1,i1) = sqrt(s);
           else
               sum2(i2,i1) = s/sum2(i1,i1);
           end
       end
   end
   
   for i1 = 1:p
       for i2 = 1:i1
           if i1 == i2
               s = 1;
           else
               s = 0;
           end
           for i3 = i1-1:-1:i2
               s = s-sum2(i1,i3)*ichol(i2,i3);
           end
           ichol(i2,i1) = s/sum2(i1,i1);
       end
   end
   
   for i1 = p:-1:1
       for i2 = 1:i1
           if i1<i2
               s = 0;
           else
               s = ichol(i2,i1);
           end
           for i3 = i1+1:p
               s = s-sum2(i3,i1)*ichol(i2,i3);
           end
           ichol(i1,i2) = s/sum2(i1,i1);
           ichol(i2,i1) = s/sum2(i1,i1);
       end
   end
   
   for i1 = 1:p
       for i2 = 1:nRow
           s = sum(sum1(i2,:)*ichol(:,i1));
           sgn = sign(s);
           if sgn>0 || non_negative<=0
               t = abs(psi(i2)*alpha*((epsl+abs(s))^spl));
               if abs(s)>t
                   L(i2,i1) = s-sgn*t;
               else
                   L(i2,i1) = 0;
               end
           else
               L(i2,i1) = 0;
           end
       end
   end
   
   if nL>0 && nL<p && cyc>bL
       for i1 = 1:nRow
           e_ssxx_n(1:nL) = -1;
           for i2 = 1:p
               s = abs(L(i1,i2));
               if s>e_ssxx_n(nL)
                   i3 = nL;
                   while i3>1 && s>e_ssxx_n(i3-1)
                       e_ssxx_n(i3) = e_ssxx_n(i3-1);
                       i3 = i3-1;
                   end
                   e_ssxx_n(i3) = s;
               end
           end
           s = e_ssxx_n(nL);
           for i2 = 1:p
               if s>abs(L(i1,i2))
                   L(i1,i2) = 0;
               end
           end
       end
   end
   
   if lL>0 && lL<nRow && cyc>bL
       for i2 = 1:p
           psi(1:lL) = -1;
           for i1 = 1:nRow
               s = abs(L(i1,i2));
               if s>psi(lL)
                   i3 = lL;
                   while i3>1 && s>psi(i3-1)
                       psi(i3) = psi(i3-1);
                       i3 = i3-1;
                   end
                   psi(i3) = s;
               end
           end
           s = psi(lL);
           for i1 = 1:nRow
               if s>abs(L(i1,i2))
                   L(i1,i2) = 0;
               end
           end
       end
   end
   
   t = 0;
   for i1 = 1:nRow
       s = sum(L(i1,:).*sum1(i1,:));
       if abs(s)>t
           t = abs(s);
       end
       psi(i1) = xx(i1) - in*s;
       if psi(i1)<eps
           psi(i1) = eps;
       end
   end
   
   if t<eps
       psi(:) = eps;
       lapla(:,:) = eps;
       break;
   end
   
   if scale>0
       for i2 = 1:p
           s = sum(L(:,i2).*L(:,i2));
           s = s*in;
           s = sqrt(s)+epsl;
           s = scale/s;
           L(:,i2) = L(:,i2)*s;
           s = s*s;
           s = s^spz;
           lapla(:,i2) = lapla(:,i2)*s;
       end
   end
   
end

%disp(cycle);
Z = zeros(p,nCol);
if t>=eps
    for i1 = 1:p
        Lpsi(:,i1) = L(:,i1)./psi;
    end
    for i1 = 1:p
        LpsiL(i1,i1) = sum(Lpsi(:,i1).*L(:,i1));
    end
    for i1 = 1:p-1
        for i3 = i1+1:p
            LpsiL(i1,i3) = sum(Lpsi(:,i1).*L(:,i3));
            LpsiL(i3,i1) = sum(Lpsi(:,i1).*L(:,i3));
        end
    end
    
    for j = 1:nCol
        tLpsiL = LpsiL;
        
        for i1 = 1:p
            tLpsiL(i1,i1) = tLpsiL(i1,i1)+lapla(j,i1);
        end
        
        for i1 = 1:p
            for i2 = i1:p
                s = tLpsiL(i1,i2);
                for i3 = i1-1:-1:1
                    s = s-tLpsiL(i1,i3)*tLpsiL(i2,i3);
                end
                ichol(i1,i2) = 0;
                ichol(i2,i1) = 0;
                if i1==i2
                    tLpsiL(i1,i1) = sqrt(s);
                else
                    tLpsiL(i2,i1) = s/tLpsiL(i1,i1);
                end
            end
        end
        
        for i1 = 1:p
            for i2 = 1:i1
                if i1 == i2 
                    s = 1; 
                else
                    s = 0;
                end
                for i3 = i1-1:-1:i2
                    s = s-tLpsiL(i1,i3)*ichol(i2,i3);
                end
                ichol(i2,i1) = s/tLpsiL(i1,i1);
            end
        end
        
        for i1 = p:-1:1
            for i2 = 1:i1
                if i1<i2 
                    s = 0;
                else
                    s = ichol(i2,i1); 
                end
                for i3 = i1+1:p
                    s = s-tLpsiL(i3,i1)*ichol(i2,i3);
                end
                ichol(i1,i2) = s/tLpsiL(i1,i1);
                ichol(i2,i1) = s/tLpsiL(i1,i1);
            end
        end
        
        for i3 = 1:p
            e_ssxx_n(i3) = sum(Lpsi(:,i3).*matrix(:,j));
        end
        
        for i1 = 1:p
            t = sum(ichol(i1,:).*e_ssxx_n);
            if t<epsl && non_negative>0
                t = 0;
            end
            Z(i1,j) = t;
        end
    end
end

resfabic = struct('Z',Z,'L',L,'psi',psi,'lapla',lapla);
resinfo = infocontent(resfabic,matrix,nRow,nCol,p);
result = extractbic(resinfo,matrix,nRow,nCol,p);
end
