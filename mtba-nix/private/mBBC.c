/*-------------------------------------------------------------------------------------
 * Modified version of func_BBC.c by JiaJun Gu
 * Changes required to make the program work under matlab environment were done.
 * 2013 Sumanik Singh <sumaniksingh@gmail.com>
 *-------------------------------------------------------------------------------------*/

/* This file contains various functions for Bayesian BiClustering.

 * The Bayesian BiClustering (BBC) program was written by Jiajun Gu
 * at the School of Engineering and Applied Sciences at Harvard University.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 *
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 *
 */
#include <stdio.h>

#include <mex.h>
#include "hBBC.h"

void mexFunction(int nlhs,mxArray *plhs[],
                     int nrhs, const mxArray *prhs[])
{
 float **tau_a, **tau_b,*tau_e;
 int ***delta, ***kapa;
 int **sdelta,**skapa;
 double *rowxNum,*colxNum;
 char **ym;
 float **ek, **bk, *qq, *smu,**salpha,**sbeta, *taua_a, *taua_b, *taub_a, *taub_b, *ek_a, *ek_b, *bk_a, *bk_b, taue_a, taue_b ; 
 double *y,*temp;
 float **matrix;
 float bic[2],normrange;
 int n,m,nk,normchoice, normr,iter,n1,n2;
 y=mxGetPr(prhs[0]);
 n=mxGetM(prhs[0]);
 m=mxGetN(prhs[0]);
 nk = (int) mxGetScalar(prhs[1]);
 normchoice = (int) mxGetScalar(prhs[2]);
 normr = mxGetScalar(prhs[3]);
 normrange = (float)normr/(float)100;
 if(normr>100 || normr<=5 && (normchoice==3 || normchoice==4))
                  mexErrMsgTxt("Wrong value for alpha parameter for normalisation. Alpha should be between 5 and 100.");
                  
 tau_a=mxCalloc(Nss,sizeof(float*));
  tau_b=mxCalloc(Nss,sizeof(float*));
  bk=mxCalloc(Nss,sizeof(float*));
  ek=mxCalloc(Nss,sizeof(float*));
  tau_e=mxCalloc(Nss,sizeof(float));
  
  qq=mxCalloc(nk,sizeof(float));
  delta=mxCalloc(Nss,sizeof(int**));
  kapa=mxCalloc(Nss,sizeof(int**));
  sdelta=mxCalloc(nk,sizeof(int*));
  skapa=mxCalloc(nk,sizeof(int*));
  smu=mxCalloc(nk,sizeof(float));
  salpha=mxCalloc(nk,sizeof(float));
  sbeta=mxCalloc(nk,sizeof(float));
  taua_a=mxCalloc(nk,sizeof(float));
  taua_b=mxCalloc(nk,sizeof(float));
  taub_a=mxCalloc(nk,sizeof(float));
  taub_b=mxCalloc(nk,sizeof(float));
  bk_a=mxCalloc(nk,sizeof(float));
  bk_b=mxCalloc(nk,sizeof(float));
  ek_a=mxCalloc(nk,sizeof(float));
  ek_b=mxCalloc(nk,sizeof(float));
  
  for (int i=0;i<nk;i++){
    salpha[i]=mxCalloc(n,sizeof(float));
    sbeta[i]=mxCalloc(m,sizeof(float));
    sdelta[i]=mxCalloc(n,sizeof(int));
    skapa[i]=mxCalloc(m,sizeof(int));
    taua_a[i]=10;
    taua_b[i]=1;
    taub_a[i]=10;
    taub_b[i]=1;
    bk_a[i]=10;
    bk_b[i]=1;
    ek_a[i]=10;
    ek_b[i]=1;
  }
  taue_a=10;
  taue_b=1;
  for(int i=0;i<Nss;i++){
    bk[i]=mxCalloc(nk,sizeof(float));
    ek[i]=mxCalloc(nk,sizeof(float));
    delta[i]=mxCalloc(nk,sizeof(int*));
    kapa[i]=mxCalloc(nk,sizeof(int*));
    tau_a[i]=mxCalloc(nk,sizeof(float));
    tau_b[i]=mxCalloc(nk,sizeof(float));
    for(int ik=0;ik<nk;ik++){
      delta[i][ik]=mxCalloc(n,sizeof(int));
      kapa[i][ik]=mxCalloc(m,sizeof(int));
    }
  }
 matrix = mxCalloc(n,sizeof(float*));
 
 for(int i=0;i<n;i++){
         temp = y;
         matrix[i] = mxCalloc(m,sizeof(float*));
         for(int j=0;j<m;j++){
                 matrix[i][j]=*y;     
                 y=y+n;
         } 
         y = temp+1 ;    
 }
 
 switch (normchoice){
  case 1:
    rowstd(n,m,matrix);//row normalization
    printf("*** Normalization finished ***\n");
    break;
  case 2:
    colstd(n,m,matrix);//column normalization
    printf("*** Normalization finished ***\n");    
    break;
  case 3:
    norm_iqr(n,m,matrix,&normrange);//iqr
    printf("*** Normalization finished ***\n");    
    break;
  case 4:
    norm_sqr(n,m,matrix,&normrange);//sqrn
    printf("*** Normalization finished ***\n");
    break;
  default:
    break;
  }
 for(iter=0;iter<nk;iter++)
     qq[iter]=0.1;
 plhs[0]=mxCreateDoubleMatrix(nk,n,mxREAL);
 plhs[1]=mxCreateDoubleMatrix(nk,m,mxREAL);
 rowxNum=(double*)mxGetPr(plhs[0]);
 colxNum=(double*)mxGetPr(plhs[1]);
 gibbssample(nk,n, m,matrix,bk,tau_a,tau_b, tau_e,ek,
	      bk_a, bk_b, taua_a,taua_b,taub_a,taub_b, ek_a, ek_b, taue_a, taue_b,
	      qq,delta,kapa,sdelta,skapa,smu,salpha,sbeta,bic); 
 
 int ik,i,j,j1;
 int *isempty;
  isempty=mxCalloc(nk,sizeof(int));
 
    mxFree(tau_a); mxFree(tau_b); mxFree(bk); mxFree(ek); mxFree(tau_e); mxFree(qq); mxFree(delta);
    mxFree(kapa); mxFree(isempty); mxFree(smu); mxFree(salpha); mxFree(sbeta); mxFree(taua_a);
    mxFree(taua_b); mxFree(taub_a); mxFree(taub_b); mxFree(bk_a); mxFree(bk_b); mxFree(ek_a); mxFree(ek_b);
 for (ik=0;ik<nk;ik++){
    n1=0;n2=0;isempty[ik]=1;
    for (i=0;i<n;i++){
      if (sdelta[ik][i]>0) n1++;
    }
    for (j=0;j<m;j++){
      if(skapa[ik][j]>0) n2++;
    }
    if (n1>1&&n2>1){
      j1++;
      isempty[ik]=0;
    }

  }
   double *temp1,*temp2;
   for(ik=0;ik<nk;ik++){
       temp1=rowxNum;
       temp2=colxNum;
     for(i=0;i<n;i++){
         *rowxNum = sdelta[ik][i];
         rowxNum=rowxNum+nk;
     }
       rowxNum=temp1+1;
     for(j=0;j<m;j++){
         *colxNum = skapa[ik][j];
         colxNum=colxNum+nk;
     }
       colxNum=temp2+1;
}
   
//    for(int i=0;i<n;i++){
//          temp = y;
//          matrix[i] = mxCalloc(m,sizeof(float*));
//          for(int j=0;j<m;j++){
//                  matrix[i][j]=*y;     
//                  y=y+n;
//          } 
//          y = temp+1 ;    
//  }
}
