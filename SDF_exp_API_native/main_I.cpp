#include <SDF_I.h>

int main(int argc, char **argv)
{
  int N=8;

  double* vecd=new double[N]();
  //  char* file_name = "../SDF_Arnoldi_native/elements_5_orbitals.txt";

  char* file_name = "../hamiltonian.txt";

  SDF_SPMV_init_txt(N, vecd, file_name);
}
