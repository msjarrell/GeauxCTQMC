#include "SDF.h"

int main()
{
  //----------------------------prologue-----------------------------//
  int N=1024;

  double* A=new double[N*N]();

  // vector to store main diagonal
  double* vecd=new double[N]();

  ///loading the matrix///

  char * pEnd;
  int row;
  int col;
  double val;

  char mystring[100];
  const char* const file_name = "elements_5_orbitals.txt";

  FILE *f = fopen(file_name, "r");

  if ( !f )
    {
      fprintf(stderr,"Could not open %s for output (%s).\n",
              file_name, strerror(errno));
      exit(1);
    }

  while( fgets(mystring, sizeof(mystring), f) ){
    val = strtod(mystring, &pEnd);
    row = strtol(pEnd, &pEnd, 10);
    col = strtol(pEnd, &pEnd, 10);

    A[row+col*N]=val;

    if(row==col)
      vecd[row]=val;
  }

  //----------------------------end of prologue-----------------------------//

  int m=5;
  double Ctemp=-0.1;
  double* v2=new double[N]();
  double t0, t1;
  double* inVec=new double[N]();

  for (int i=0; i<N; i++){
    inVec[i]=rand();
  }

  t0=time_wall_fp();
  {
    mat_exp_krylov(N, m, Ctemp, inVec, vecd, v2);
  }
  t1=time_wall_fp();

  printf("total elapsed time using SDF_krylov is: %.6f s\n",
         (t1-t0));
  
  double* v1=new double[N]();
  double* expA=new double[N*N]();

  t0=time_wall_fp();
  {
    mat_exp(A, expA, N, Ctemp);
    matmul(expA, inVec, v1, N, 1, N);
  }
  t1=time_wall_fp();
  printf("total elapsed time with direct approach is: %.3f s\n",
         (t1-t0));

  double ERR=0.0;
  for(int i=0; i<N; i++)
    {
      ERR+= fabs( v1[i]-v2[i] );
    }

  printf("   ERR = %f \n", ERR);

  double sos = 0;
  for (int i=0; i<N; i++)
    {
      double d = inVec[i]=rand();
      sos += d * d;
    }

  double nf = pow(sos,-0.5);
  for (int i=0; i<N; i++) inVec[i] *= nf;

  t0=time_wall_fp();
  {
    mat_exp(A, expA, N, Ctemp);
    matmul(expA, inVec, v1, N, 1, N);
  }
  t1=time_wall_fp();
  printf("total elapsed time with direct approach is: %.3f s\n",
         (t1-t0));

  ERR=0.0;
  for(int i=0; i<N; i++)
    {
      ERR+= fabs( v1[i]-v2[i] );
    }

  printf("   ERR = %f \n", ERR);

}
