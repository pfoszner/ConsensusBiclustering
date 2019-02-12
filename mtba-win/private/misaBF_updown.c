/*-----------------------------------------------------------------------------
 * misaBF_updown.c:
 * Author: Jayesh Kumar Gupta, 2013.
 * Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
 *          Indian Institute of Technology, Kanpur, India
 *
 *---------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include "mex.h"

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
  double *data, mean, sd, absmax, eff_thr_up, eff_thr_down, tmp;
  double thr;

  plhs[0] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  

  data = mxGetPr(prhs[0]);
  
  for (i=0; i<ncol; i++){
    thr = threshold[i];

    /* column dependent thresholds */

    mean=0.0; sd=0.0; absmax=0.0;

    for (j=0; j<nrow; j++){
      sd += j*pow(*data-mean, 2)/(j+1);
      mean += (*data-mean)/(j+1);
      ++data;
    }
    sd = sqrt(sd/(nrow-1));
    eff_thr_up = mean + thr*sd;
    eff_thr_down = mean - thr*sd;

    /* threshold + search for maximum */

    data -= nrow;
    for (j=0; j<nrow; j++){
      if (*data <= eff_thr_up && *data >= eff_thr_down){
        *data = 0;
      } else {
        tmp = fabs(*data);
        if (tmp > absmax) { absmax = tmp; }
      }
      data++;
    }

    /* norm(1) */

    if (absmax != 0){
      data -= nrow;
      for (j=0; j<nrow; j++){
        *data /= absmax;
        data++;
      }
    }
  }
  /* mxSetPr(plhs[0], result); */
}
