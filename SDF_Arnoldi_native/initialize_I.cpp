#include "SDF_I.h"

using namespace std;

int main()
{
  int N=1024;

  double* vecd=new double[N]();

  char* file_name = "../SDF_Arnoldi_native/elements_5_orbitals.txt";

  SDF_SPMV_init_txt(N, vecd, file_name);
  
}
