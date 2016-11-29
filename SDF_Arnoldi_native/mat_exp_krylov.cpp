#include "SDF.h"

using namespace std;

void mat_exp_krylov(int N, int m, double Ctemp, double* inVec, double* vecd, double* outVec)
//                    outVec = exp ( Ctemp * norm(inVec) * SPM )
// N is the dimension of SPM, m is the desired dimension of the krylov subspace
// outVec is an empty array of size N, vecd is the main diagonal of SPM
// inVec is the input array, this function will normalize inVec before applying it.
{
#ifdef TIMING
  double t_start = time_wall_fp();
#endif
  double* H=new double[m*m]();
  double* Vmat=new double[N*m];
  double* expH=new double[m*m]();

  double normal;
  normal=norm(inVec, N);
  for (int i=0; i<N; i++){
    inVec[i]=inVec[i]/normal;
  }
#ifdef TIMING
  double t_prep = time_wall_fp();
#endif
  arnoldi_SDF(Vmat, inVec, vecd, m, H, N);
#ifdef TIMING
  double t_arnoldi = time_wall_fp();
#endif
  mat_exp(H, expH, m, Ctemp);
#ifdef TIMING
  double t_exp = time_wall_fp();
#endif
  matmul(Vmat, expH, outVec, N, 1, m);
#ifdef TIMING
  double t_end = time_wall_fp();
  printf("Timings: %.3f ms = %.3f + %.3f + %.3f + %.3f  (prep,arn,exp,mm) \n",
         ( t_end - t_start ) * 1e3,
         ( t_prep - t_start ) * 1e3,
         ( t_arnoldi - t_prep ) * 1e3,
         ( t_exp - t_arnoldi ) * 1e3,
         ( t_end - t_exp ) * 1e3);
#endif

  delete H;
  delete Vmat;
  delete expH;
}
