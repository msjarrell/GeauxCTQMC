#ifndef SDF_D
#define SDF_D



#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>
#include <errno.h>
#include <cmath>




#ifdef HAVE_MKL

#include <mkl.h>

#else

extern "C" void dsyev_
( const char *jobz, const char *uplo, int *n, double *a, int *lda, 
  double *w, double *work, int *lwork, int *info );

#define dsyev dsyev_

#endif

using namespace std;

void arnoldi_SDF(double* Vmat, double* inVec, double* vecd, int m, double* H, int N);

extern "C" void SDF_SPMV_p(double *C, double *B);
extern "C" void SDF_SPMV_n(double *C, double *B);

void mat_exp_krylov(int N, int m, double Ctemp, double* inVec, double* vecd, double* outVec);

inline double
time_wall_fp()
{
  struct timespec now;
  clock_gettime(CLOCK_REALTIME,&now);
  return now.tv_sec + ((double)now.tv_nsec) * 1e-9;
}

void mat_exp(double* Ain, double* B, int N, double t);

double norm(double* x, int N);

double dot(double* x, double *y, int N);

void matmul(double* x, double* y, double* z, int m, int n, int p);

#endif
