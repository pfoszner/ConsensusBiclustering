/*-----------------------------------------------------------------------------
 * misaBF_up.c:
 * Author: Jayesh Kumar Gupta, 2013.
 * Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
 *          Indian Institute of Technology, Kanpur, India
 *
 *---------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include "mex.h"
/* #include "matrix.h" */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
  if (nrhs!=2){
    mexErrMsgTxt("Wrong number of inputs");
  }

  mwSize nrow = mxGetM(prhs[0]);
  mwSize ncol = mxGetN(prhs[0]);
  double *threshold;
  threshold = mxGetPr(prhs[1]);

  mwSize i, j;
  double *data, mean, sd, absmax, eff_thr, eff_max;
  double thr;

  plhs[0] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);

  double *result;

  data = mxGetPr(prhs[0]);
  result = mxGetPr(plhs[0]);

  for (i=0; i<ncol; i++){

    thr = threshold[i];

    /* column dependent thresholds, plus search for absolute maximum */

    mean=0.0; sd=0.0; absmax=0.0;

    for (j=0; j<nrow; j++){
      sd += j*pow(*data-mean, 2)/(j+1);
      mean += (*data-mean)/(j+1);
      if(*data > absmax) { absmax = *data; }
      ++data;
    }
    sd = sqrt(sd/(nrow-1));
    eff_max = mean + thr*sd;
    eff_max = absmax;

    /* threshold + norm (1)*/
    data -= nrow;

    for (j=0; j<nrow; j++){
      if (*data <= eff_thr){
        *data = 0;
      } else{
        *data /= eff_max;
      }
      data++;
    }
  }
  /* mxSetData(plhs[0], result); */
  memcpy(result, data, sizeof(double)*nrow*ncol);
}
