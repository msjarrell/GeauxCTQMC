#include "SDF_I.h"

using namespace std;

void SDF_SPMV_init(int* row, int* col, double* val, int N, int NNZ, double* vecd , double J)
// inputs should be the sparse matrix in COO format, sorted on a row major format (assuming symmetry, row/column major are the same!)
// N is the dimension of matrix and NNZ is the number of nonzero elements
// vecd should be an ampty array of size N, it will return the main diagonal
{
  //----------------------------prologue-----------------------------//

  // CSR to store positive Js
  double* CSR_D_p = new double[NNZ]();
  int* CSR_Cid_p = new int[NNZ]();
  int* CSR_Rpt_p = new int[N+1]();

  // CSR to store negative Js
  double* CSR_D_n = new double[NNZ]();
  int* CSR_Cid_n = new int[NNZ]();
  int* CSR_Rpt_n = new int[N+1]();

  int k_p=0;
  int l_p=0;
  int k_n=0;
  int l_n=0;

  CSR_Rpt_p[0]=0;
  l_p=1;
  CSR_Rpt_n[0]=0;
  l_n=1;

  for (int ii=0; ii<NNZ; ii++){
    if(row[ii]==col[ii]){
      vecd[row[ii]]=val[ii];
    } else {
      if(val[ii]>0){
	CSR_D_p[k_p] = val[ii];
	CSR_Cid_p[k_p] = col[ii];
	while(l_p<=row[ii]){
	  CSR_Rpt_p[l_p]=k_p;
	  l_p++;
	}
	k_p++;
      } else {
	CSR_D_n[k_n] = val[ii];
	CSR_Cid_n[k_n] = col[ii];
	while(l_n<=row[ii]){
	  CSR_Rpt_n[l_n]=k_n;
	  l_n++;
	}
	k_n++;
      }
    }
  }
  for (int i=l_n; l_n<N+1; i++){
    CSR_Rpt_n[l_n]=k_n;
    l_n++;
  }
  for (int i=l_p; l_p<N+1; i++){
    CSR_Rpt_p[l_p]=k_p;
    l_p++;
  }

  double **DIA_D_p = new double*[N]();
  int* DIA_Did_p = new int[N]();
  double **DIA_D_n = new double*[N]();
  int* DIA_Did_n = new int[N]();

  CSR2DIA(CSR_D_p, CSR_Cid_p, CSR_Rpt_p, CSR_Rpt_p[N], DIA_D_p, DIA_Did_p, N);

  CSR2DIA(CSR_D_n, CSR_Cid_n, CSR_Rpt_n, CSR_Rpt_n[N], DIA_D_n, DIA_Did_n, N);


  ///DIA to SDF///

  int NNZD_p=0;
  for(int i=0; i<N; i++){
    if(DIA_Did_p[i]!=0){
      NNZD_p++;
    }
  }

  int **SDF_p = new int*[NNZD_p];
  DIA2SDF(DIA_D_p, DIA_Did_p, N, SDF_p);

  int NNZD_n=0;
  for(int i=0; i<N; i++){
    if(DIA_Did_n[i]!=0){
      NNZD_n++;
    }
  }

  int **SDF_n = new int*[NNZD_n];
  DIA2SDF(DIA_D_n, DIA_Did_n, N, SDF_n);


  ///SDF_SPMV initialization///

  dyn_compile_p(SDF_p, NNZD_p, N, J);
  dyn_compile_n(SDF_n, NNZD_n, N, J);

  delete CSR_D_p;
  delete CSR_Cid_p;
  delete CSR_Rpt_p;
  delete CSR_D_n;
  delete CSR_Cid_n;
  delete CSR_Rpt_n;
  delete DIA_D_p;
  delete DIA_Did_p;
  delete DIA_D_n;
  delete DIA_Did_n;
  delete SDF_p;
  delete SDF_n;

}
