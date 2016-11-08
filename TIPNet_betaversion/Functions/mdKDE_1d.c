#include <stdio.h>
#include <math.h>
#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  double *X; 
  double *minind, *maxind; 
  double *ICoords, *N, *h;
  double *pdf;
  int  d, ii,ic,nI,nval,minI,maxI, n[3],nn;
  double prodi, ui;

  X         = (double *)mxGetPr(prhs[0]); //pointer to first input: data array
  d         = mxGetN(prhs[0]); // getM: number of rows
  minind    = (int *)mxGetPr(prhs[1]); //pointer to second input: min indices
  maxind    = (int *)mxGetPr(prhs[2]); //pointer to third input: max indices
  ICoords   = (double*)mxGetPr(prhs[3]); //pointer to 4th input: i coordinates
        n[0]=mxGetN(prhs[3]);
  N =        (int *)mxGetPr(prhs[4]); //point to total number of coords for pdf
  h =        (double *)mxGetPr(prhs[5]);
  
  //printf("n values = %d  %d  %d\n",n[0],n[1],n[2]);

  nval = *N;
  
  plhs[0] = mxCreateNumericMatrix(1,nval, mxDOUBLE_CLASS, mxREAL);
  pdf = mxGetPr(plhs[0]);
  
  // printf("N = %d\n",nval);
  // printf("dim = %d\n",d);
  
  for (ii=0; ii<n[0]; ii++){
      //printf("iteration = %d \n",ii);
      ui = (X[0]-ICoords[ii]) / h[0];
      if (ui*ui < 1){
          ic = ii + minind[0] - 1;
          prodi = 0.75 * (1-(ui*ui));
          pdf[ic] = prodi;}
  }

//       free(ICoords); 
//       free(JCoords);
//       free(KCoords);
//       free(N);
   mxSetPr(plhs[0], pdf);
  // memmove(mxGetData(plhs[0]),pdf,nval*nval*nval*sizeof(double)); //memmove(destination, source, size)
   // memmove(plhs[0],pdf,nval*nval*nval*sizeof(double));
   //mxFree(pdf);

//    free(N);


}
