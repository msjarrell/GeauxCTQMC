#include <SDF.h>

int main(int argc, char **argv)
{
  //----------------------------prologue-----------------------------//
  int N=1024;

  // vector to store main diagonal
  double* vecd=new double[N]();

  char * pEnd;
  double val;

  char mystring[100];
  const char* const file_name = "vecd.txt";

  FILE *f = fopen(file_name, "r");

  if ( !f )
    {
      fprintf(stderr,"Could not open %s for output (%s).\n",
              file_name, strerror(errno));
      exit(1);
    }

  int k=0;
  while( fgets(mystring, sizeof(mystring), f) ){
    vecd[k] = strtod(mystring, &pEnd);
    k++;
  }

  fclose(f);
  //----------------------------end of prologue-----------------------------//

  int m=5;
  double Ctemp=-0.1;
  double* v2=new double[N]();
  double* inVec=new double[N]();

  for (int i=0; i<N; i++){
    inVec[i]=rand();
  }

  mat_exp_krylov(N, m, Ctemp, inVec, vecd, v2);

  printf("Done!\n");

}
