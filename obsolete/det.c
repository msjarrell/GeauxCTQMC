#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <lapacke.h>

#include <complex.h>

#include "det.h"

extern void zgetrf_ (int *, int *, double complex *, int *, int *, int *);
  //void zgetri_(int* N, double complex* A, int* lda, int* IPIV, double complex* WORK, int* lwork, int* INFO);




void
det (int SZ, int N, double complex * m, double complex * determinant)
{

  double complex *result = malloc (SZ * SZ * sizeof (double complex));

  memcpy (result, m, SZ * SZ * sizeof (double complex));

  int *pivot = malloc (N * sizeof (int));
  int info;

  int i, j;

  *determinant = 1.0;

  zgetrf_ (&N, &N, result, &SZ, pivot, &info);

  for (i = 0; i < N; i++)
    (*determinant) *= result[i * SZ + i];

  free (result);
  //  zgetrf_( &N, result, &SZ, pivot,&info );  

}
