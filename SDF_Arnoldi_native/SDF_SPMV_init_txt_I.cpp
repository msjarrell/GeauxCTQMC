#include "SDF_I.h"

using namespace std;

void SDF_SPMV_init_txt(int N, double* vecd, char* file_name)
// the text file should have the sparse matrix in COO format, sorted on a row major format (assuming symmetry, row/column major are the same!)
// N is the dimension of matrix and vecd should be an ampty array of size N, it will return the main diagonal
// this function will also write out vecd to a file vecd.txt
{
  // this is our initial guesstimate for NNZ, to be safe you might want to choose larger values depending on N
  int NNZ=10000;

  int Digs = DBL_DIG;

  ///loading the matrix///

  char * pEnd;
  int* row = new int[NNZ]();
  int* col = new int[NNZ]();
  double* val = new double[NNZ]();

  double J;
  int Jset=0;

  char mystring[100];

  FILE *f = fopen(file_name, "r");

  if ( !f )
    {
      fprintf(stderr,"Could not open %s for output (%s).\n",
              file_name, strerror(errno));
      exit(1);
    }

  int ii=0;
  while( fgets(mystring, sizeof(mystring), f) ){
    val[ii] = strtod(mystring, &pEnd);
    row[ii] = strtol(pEnd, &pEnd, 10);
    col[ii] = strtol(pEnd, &pEnd, 10);

    if (row[ii]!=col[ii] & !Jset & val[ii]>0) {
      J=val[ii];
      Jset=1;
    }
    ii++;
  }
  NNZ=ii;

  fclose(f);

  SDF_SPMV_init( row, col, val, N, NNZ, vecd , J);

  const char* const output_name = "vecd.txt";

  FILE *f1 = fopen(output_name, "w");

  if ( !f1 )
    {
      fprintf(stderr,"Could not open %s for output (%s).\n",
              output_name, strerror(errno));
      exit(1);
    }

  for(int i=0; i<N; i++)
    {
      fprintf(f1, "%.*e \n", Digs, vecd[i]);
    }
  fclose(f1);

  delete row;
  delete col;
  delete val;

}
