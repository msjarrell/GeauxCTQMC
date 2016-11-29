#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>
#include <errno.h>
#include <cmath>
#include <cfloat>

void CSR2DIA(double *CSR_D, int *CSR_Cid, int *CSR_Rpt, int NNZ,
	     double **DIA_D, int *DIA_Did, int N);

void DIA2SDF(double **DIA_D, int *DIA_Did, int N, int **SDF);

void SDF_SPMV_init(int* row, int* col, double* val, int N, int NNZ, double* vecd, double J);

void SDF_SPMV_init_txt(int N, double* vecd, char* file_name);

void dyn_compile_p(int **SDF, int NNZD, int N, double J);
void dyn_compile_n(int **SDF, int NNZD, int N, double J);
