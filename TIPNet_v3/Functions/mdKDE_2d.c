#include <stdio.h>
#include <math.h>
#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  double *X; 
  double *minind, *maxind; 
  double *ICoords, *JCoords, *N, *h;
  double *pdf;
  int  d, ii,jj,kk,ic,jc,nI,nJ,nval,minI,maxI, n[3],nn;
  double prodi,prodj,prodk,prod, ui,uj,uk;

  X         = (double *)mxGetPr(prhs[0]); //pointer to first input: data array
  d         = mxGetN(prhs[0]); // getM: number of rows
  minind    = (int *)mxGetPr(prhs[1]); //pointer to second input: min indices
  maxind    = (int *)mxGetPr(prhs[2]); //pointer to third input: max indices
  ICoords   = (double*)mxGetPr(prhs[3]); //pointer to 4th input: i coordinates
        n[0]=mxGetN(prhs[3]);
  JCoords   = (double*)mxGetPr(prhs[4]); //pointer to 5th input: j coordinates
        n[1]=mxGetN(prhs[4]);

  N =        (int *)mxGetPr(prhs[5]); //point to total number of coords for pdf
  h =        (double *)mxGetPr(prhs[6]);
  
  //printf("n values = %d  %d  %d\n",n[0],n[1],n[2]);

  nval = *N;
  
  plhs[0] = mxCreateNumericMatrix(1,nval*nval, mxDOUBLE_CLASS, mxREAL);
  pdf = mxGetPr(plhs[0]);
 
 // printf("N = %d\n",nval);
 // printf("dim = %d\n",d);
  
  for (ii=0; ii<n[0]; ii++){
      //printf("iteration = %d \n",ii);
      ui = (X[0]-ICoords[ii]) / h[0];
      if (ui*ui < 1){
          for (jj=0; jj<n[1]; jj++){
              uj = (X[1]-JCoords[jj]) / h[1];
              prodi = ui*ui + uj*uj;
              if (prodi < 1){
                  ic = ii + minind[0] - 1;
                  jc = jj + minind[1] - 1;
                  pdf[jc*nval+ic] = .637 * (1-prodi);}
       
          }
    
      }
  }

   mxSetPr(plhs[0], pdf);
  // memmove(mxGetData(plhs[0]),pdf,nval*nval*nval*sizeof(double)); //memmove(destination, source, size)
   // memmove(plhs[0],pdf,nval*nval*nval*sizeof(double));
   //mxFree(pdf);

//    free(N);


}
