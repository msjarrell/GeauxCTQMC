#include <cstdio>
#include <iostream>
#include <new>

void CSR2DIA(double *CSR_D, int *CSR_Cid, int *CSR_Rpt, int NNZ,
	    double **DIA_D, int *DIA_Did, int N)
// DIA_D is an array of pointers to arrays storing the values. DIA_id is a 
// 1024 binary showing which of the diagonals is NNZ
{

  int j=0;
  int d=0;
  int k;

  for (int i=0; i<NNZ; i++) {
    while( CSR_Rpt[j] <= i ) {
      if (CSR_Rpt[j]==CSR_Rpt[j+1]) { // k is to handle zero rows
	j++;
	k++;
      } else {
	k=0;
	j++;
      }
    }
    d=CSR_Cid[i]-(j-k-1)+N-1; // diagonal index (ranges from 0 to 2046)
    if ( d<N ) {     // matrix is symethrical, only need half
      if (DIA_Did[d]==0){
	DIA_Did[d]=1; // mark that this diagonal is NNZ
	DIA_D[d]= new double[d+1];// allocate memory for the NNZ diagonal
      }
      DIA_D[d][CSR_Cid[i]]=CSR_D[i];// store based on column index
    }
  }
  return;
}
