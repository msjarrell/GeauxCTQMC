#include <stdio.h>            /* Standard Library of Input and Output */
#include <complex.h>		/* Standart Library of Complex Numbers */
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include "parameters.h"
#include "matrix.h"

#include "calc_green_func.h"


double complex
complex_exp (double t)
{
  //Assuming t is real, return e ^ (i * t)
  return (cos (t) + I * sin (t));
}



void
init_g_tau (double complex *g_tau,double beta,double v2,double D)	//initialize the table for green function
#ifndef DMFT_LOOP
//Calculate g(tau) from scratch and store the result into a table 
{
  
  double complex g_omega[2 * N_OMEGA];
  double omega;
  int i, j;
  
  //printf("%f",D);
  for (i = 0; i < 2* N_OMEGA; i++) {
    int n = 2 * (i - N_OMEGA) + 1 ;
    omega = n * PI / beta;
    g_omega[i] =  v2 * (clog(I*omega +D)-clog(I*omega-D));

  }

#pragma simd
  for (i = 0; i < N_TAU; i++) 
    g_tau[i] = 0.0 + 0.0 * I;

  for (j = 0; j < 2 * N_OMEGA; j++) {  
    double omega =  (2 * (j - N_OMEGA) + 1) * PI / beta;
#pragma simd
    for (i = 0; i < N_TAU; i++) {
      double tau = (i+0.5) * beta / (N_TAU);
      g_tau[i] += g_omega[j] * complex_exp ( -1.0 * omega * tau);

    }
  }
 
    
#pragma simd
  for (i = 0; i < N_TAU; i++){ 
    g_tau[i] = g_tau[i] / beta;
    //printf("%lf\n",g_tau[i]);
  }
}
#else
//Read g(tau) from file and store the result into a table 
{
  FILE *fp;
  fp= fopen("GS.txt","r");
  int i;
  double tr,ti;
  for(i=0;i<N_TAU;i++){
    fscanf(fp,"%lf\t%lf\n",&tr,&ti);
    g_tau[i] = tr + I * ti;
    //    printf("%lf\n",g_tau[i]);
  }
  fclose(fp);
}
#endif


void
get_g_tau (MTYPE * g_tau, double tau, double beta,MTYPE * result)
{
  
  int sign=1;

  if (tau < 0){
    tau += beta;
    sign=-1;
  }
  
  //while (tau >= beta)
  if (tau >= beta){
    tau -= beta;
    sign=-1;
  }
  //printf("tau2 %f\n",tau);
  /*
  assert(tau>=0);
  assert(tau<=beta);
  */


  int i = floor (tau * (N_TAU)/ beta + 0.5);
  int xi = i-1 + (i==0) - (i==N_TAU);
  int xj = xi + 1;

  //MTYPE gi = g_tau[xi];
  //int i = floor (tau * (N_TAU)/ beta);
  MTYPE gi = g_tau[xi];
  MTYPE gj = g_tau[xj];
  
  //MTYPE gi=-tau/beta;

  //MTYPE gt=0.1/(exp(-1.0*beta)+1)*exp(-1.0*(beta-tau));

  //*result = sign*gi;
  *result =  sign* (gi + (gj - gi) / (beta / N_TAU) * (tau - (xi+0.5) * beta / N_TAU));
  //  printf("get_g_tau done\n");
}

