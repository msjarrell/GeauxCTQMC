#include "SDF.h"


void mat_exp(double* Ain, double* B, int N, double t)
//                   B(N,N) = exp ( t * Ain(N,N) )
{
#ifdef TIMING
  double t0=time_wall_fp();
#endif

  int n = N; 
  int lda = N;
  int info;
  int lwork;
  double wkopt;
  double* w = new double[N]();
  double* A = new double[N*N]();
  // duplicate Ain so it won't get destroyed
  for (int i=0; i<N*N; i++){
    A[i] = Ain[i];
  }
  lwork = -1;

#ifdef TIMING
  double t1=time_wall_fp();
#endif
  dsyev( "Vectors", "Upper", &n, A, &lda, w, &wkopt, &lwork, &info );
  lwork = (int)wkopt;

  double* work = new double[lwork]();
#ifdef TIMING
  double t2=time_wall_fp();
#endif
  dsyev( "Vectors", "Upper", &n, A, &lda, w, work, &lwork, &info );

#ifdef TIMING
  double t3=time_wall_fp();
#endif
  for ( int i=0; i<N; i++ ) w[i] = exp(w[i]*t);
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        B[i+j*N] = B[i+j*N] + A[i+k*N] * w[k] * A[j+k*N];
      }
    }
  }
  delete w;
  delete A;
  delete work;
#ifdef TIMING
  double t4=time_wall_fp();

  printf("Timings for mat_exp: %.3f ms = %.3f + %.3f + %.3f + %.3f  (exp_prep,dsyev_prep,dsyev,exp) \n",
         ( t4 - t0 ) * 1e3,
         ( t1 - t0 ) * 1e3,
         ( t2 - t1 ) * 1e3,
         ( t3 - t2 ) * 1e3,
         ( t4 - t3 ) * 1e3);
#endif
}

double norm(double* x, int N)
{
  double norm1=0.0;
  for (int i=0; i<N; i++)
    {
      norm1 += x[i] * x[i];
    }
  norm1 = sqrt(norm1);
  return norm1;
}


double dot(double* x, double *y, int N)
{
  double dot=0.0;
  for(int i=0; i<N; i++)
    {
      dot+=x[i]*y[i];
    }
  return dot;
}

void matmul(double* x, double* y, double* z, int m, int n, int p)
//              Z(m,n) = X(m,p) * Y(p,n)
{
  for(int i=0; i<n; i++)
    {
      for(int j=0; j<m; j++)
	{
          z[ j + m*i ] = 0;
	  for(int k=0; k<p; k++)
	    {
	      z[j+m*i]+=x[j+m*k]*y[k+p*i];
	    }
	}
    }
}
