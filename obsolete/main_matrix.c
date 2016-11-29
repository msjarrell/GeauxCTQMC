#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <complex.h>

//#include "hdf5.h"

//#include "mpiprocess.h"

#include "matrix.h"
#include "parameters.h"

int
main (int argc, char **argv)
{

  double complex g_tau[N_TAU];
  double delta = 1.0;
  double beta = 10.0;
  int i, j;

  init_g_tau (g_tau, beta);

  for (i = 0; i < 1000; i++) {
    double complex g_1 = 0;
    get_g_tau (g_tau, i * 0.001 * beta, beta, &g_1);
    double complex g_2 = 0;
    calc_g_tau (i * 0.001 * beta,beta, &g_2);
    printf ("%f\t%f\t%f\t%f\n", creal (g_1), cimag (g_1), cabs (g_1),
	    cabs (g_2));
  }

  /*
     int test_sz=3;

     double **is=malloc(sizeof(double *));
     is[0]=malloc(sizeof(double)*test_sz);
     double **ie=malloc(sizeof(double *));
     ie[0]=malloc(sizeof(double)*test_sz);

     //some test numbers
     is[0][0]=0.047235;
     is[0][1]=0.139749;
     is[0][2]=0.799764;
     ie[0][0]=0.137529;
     ie[0][1]=0.306733;
     ie[0][2]=0.833604;

     Matrix *m=malloc(sizeof(Matrix));

     construct_matrix( is, ie, &test_sz, m, g_tau, beta);

     int sz=m->max_sz;
     for ( i=0;i<test_sz;i++){
     printf("{");
     for ( j=0;j<test_sz;j++)
     printf("%f+%f*I,\t",creal(m->mat[i*sz+j]),cimag(m->mat[i*sz+j]));
     printf("},\n");
     }

     double complex result=m->det;

     printf("The determinant is %f + %f *I\n",creal(result),cimag(result));
     printf("The mode of determinant is %g\n",cabs(result));
   */
  return 0;
}
