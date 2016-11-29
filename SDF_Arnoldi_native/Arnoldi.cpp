#include "SDF.h"

using namespace std;

void arnoldi_SDF(double* Vmat, double* inVec, double* vecd, int m, double* H, int N)
{

  double* V=new double[N];

  double normal=norm(inVec, N);
  double INV;

  for(int i=0; i<N; i++)
    {
      Vmat[i] = inVec[i] / normal;
    }

  for(int j = 0; j < m-1; j++)
    {
      for (int i=0; i<N; i++)
	{
	  V[i]=vecd[i]*Vmat[j*N+i];
	}
      SDF_SPMV_p(V, &Vmat[j*N]);
      SDF_SPMV_n(V, &Vmat[j*N]);
      for(int i = 0; i <= j; i++)
	{
          H[i+j*m]=dot(V, &Vmat[i*N], N);
	  for(int k = 0; k < N; k++)
	    {
	      V[k]=V[k]-H[i+j*m]*Vmat[k+i*N];
	    }
	}
      H[j+1+j*m]=norm(V, N);
      INV=(1.0/H[j+1+j*m]);
      for(int i = 0; i < N; i++)
	{
	  Vmat[i+(j+1)*N]=V[i]*INV;
	}
    }
  int j = m-1;
  for (int i=0; i<N; i++)
    {
      V[i]=vecd[i]*Vmat[j*N+i];
    }
  SDF_SPMV_p(V, &Vmat[j*N]);
  SDF_SPMV_n(V, &Vmat[j*N]);
  for(int i = 0; i <= j; i++)
    {
      H[i+j*m]=dot(V, &Vmat[i*N], N);
      for(int k = 0; k < N; k++)
	{
	  V[k]=V[k]-H[i+j*m]*Vmat[k+i*N];
	}
    }

  delete V;
  return;
}
