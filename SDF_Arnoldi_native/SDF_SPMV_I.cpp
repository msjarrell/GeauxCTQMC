#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>
#include <errno.h>
#include <dlfcn.h>
#include <cfloat>

using namespace std;

void dyn_compile_p(int **SDF, int NNZD, int N, double J)
{
  // Name of file that will contain the SDF_SPMV routine body.
  //
  const char* const file_name = "SDF_SPMV_dc_func_p.cpp";

  int Digs = DBL_DIG;

  //
  // Write out the SDF_SPMV body.
  //

  FILE *f = fopen(file_name, "w");

  if ( !f )
    {
      fprintf(stderr,"Could not open %s for output (%s).\n",
              file_name, strerror(errno));
      exit(1);
    }

  fprintf(f," extern \"C\" void\n");
  fprintf(f," SDF_SPMV_p(double *C, double *B)\n");
  fprintf(f," //                 out_vec    in_vec\n");
  fprintf(f," {\n");
  fprintf(f,"   int current_position, k;\n");

  // each NNZ diagonal will have it's own loop body.
  for (int i0=0; i0<NNZD; i0++)
    {
      //base (first NNZ)
      fprintf(f," k=%d;\n", SDF[i0][1]);

      fprintf(f," for(int j=0; j<%d; j++){\n", N/SDF[i0][3]);

      // body of N nested loops where N=depth
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," for (int i%d=0; i%d<%d; i%d++){\n", i1, i1, SDF[i0][3+SDF[i0][2]+i1], i1);
	}

      // index of the next NNZ
      fprintf(f,"current_position=k");
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f,"+i%d*%d", i1, SDF[i0][4+i1]);
	}
      fprintf(f,";\n");

      // if not out of bound do the computation
      fprintf(f," if (current_position<%d){\n", SDF[i0][0]+N);

      fprintf(f,"C[current_position+%d] += %.*e * B[current_position+%d];\n", -SDF[i0][0], Digs, J, -SDF[i0][0]);

      fprintf(f," }\n");

      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," }\n");
	}

      // move base to next location
      fprintf(f,"k+=%d;\n", SDF[i0][3]);

      fprintf(f,"}\n");

    }

  // this version of SDF assumes symmetry, the following will handle the upper half of the matrix.
  if(SDF[NNZD-1][0]!=0){
      //base (first NNZ)
      fprintf(f," k=%d;\n", SDF[NNZD-1][1]);

      fprintf(f," for(int j=0; j<%d; j++){\n", N/SDF[NNZD-1][3]);

      // body of N nested loops where N=depth
      for (int i1=0; i1<SDF[NNZD-1][2]-1; i1++)
	{
	  fprintf(f," for (int i%d=0; i%d<%d; i%d++){\n", i1, i1, SDF[NNZD-1][3+SDF[NNZD-1][2]+i1], i1);
	}

      // index of the next NNZ
      fprintf(f,"current_position=k");
      for (int i1=0; i1<SDF[NNZD-1][2]-1; i1++)
	{
	  fprintf(f,"+i%d*%d", i1, SDF[NNZD-1][4+i1]);
	}
      fprintf(f,";\n");

      // if not out of bound do the computation
      fprintf(f," if (current_position<%d){\n", SDF[NNZD-1][0]+N);

      fprintf(f,"C[current_position] += %.*e * B[current_position];\n", Digs, J);

      fprintf(f," }\n");

      for (int i1=0; i1<SDF[NNZD-1][2]-1; i1++)
	{
	  fprintf(f," }\n");
	}

      // move base to next location
      fprintf(f,"k+=%d;\n", SDF[NNZD-1][3]);

      fprintf(f,"}\n");

    }

  for (int i0=0; i0<NNZD-1; i0++)
    {
      //base (first NNZ)
      fprintf(f," k=%d;\n", SDF[i0][1]);

      fprintf(f," for(int j=0; j<%d; j++){\n", N/SDF[i0][3]);

      // body of N nested loops where N=depth
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," for (int i%d=0; i%d<%d; i%d++){\n", i1, i1, SDF[i0][3+SDF[i0][2]+i1], i1);
	}

      // index of the next NNZ
      fprintf(f,"current_position=k");
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f,"+i%d*%d", i1, SDF[i0][4+i1]);
	}
      fprintf(f,";\n");

      // if not out of bound do the computation
      fprintf(f," if (current_position<%d){\n", SDF[i0][0]+N);

      fprintf(f,"C[current_position] += %.*e * B[current_position];\n", Digs, J);

      fprintf(f," }\n");

      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," }\n");
	}

      // move base to next location
      fprintf(f,"k+=%d;\n", SDF[i0][3]);

      fprintf(f,"}\n");

    }

  fprintf(f," return;\n");
  fprintf(f," }\n");

  fclose(f);

}


void dyn_compile_n(int **SDF, int NNZD, int N, double J)
{
  // Name of file that will contain the SDF_SPMV routine body.
  //
  const char* const file_name = "SDF_SPMV_dc_func_n.cpp";

  int Digs = DBL_DIG;

  //
  // Write out the SDF_SPMV body.
  //

  FILE *f = fopen(file_name, "w");

  if ( !f )
    {
      fprintf(stderr,"Could not open %s for output (%s).\n",
              file_name, strerror(errno));
      exit(1);
    }

  fprintf(f," extern \"C\" void\n");
  fprintf(f," SDF_SPMV_n(double *C, double *B)\n");
  fprintf(f," //                 out_vec    in_vec\n");
  fprintf(f," {\n");
  fprintf(f,"   int current_position, k;\n");

  // each NNZ diagonal will have it's own loop body.
  for (int i0=0; i0<NNZD; i0++)
    {
      //base (first NNZ)
      fprintf(f," k=%d;\n", SDF[i0][1]);

      fprintf(f," for(int j=0; j<%d; j++){\n", N/SDF[i0][3]);

      // body of N nested loops where N=depth
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," for (int i%d=0; i%d<%d; i%d++){\n", i1, i1, SDF[i0][3+SDF[i0][2]+i1], i1);
	}

      // index of the next NNZ
      fprintf(f,"current_position=k");
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f,"+i%d*%d", i1, SDF[i0][4+i1]);
	}
      fprintf(f,";\n");

      // if not out of bound do the computation
      fprintf(f," if (current_position<%d){\n", SDF[i0][0]+N);

      fprintf(f,"C[current_position+%d] += -%.*e * B[current_position+%d];\n", -SDF[i0][0], Digs, J, -SDF[i0][0]);

      fprintf(f," }\n");

      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," }\n");
	}

      // move base to next location
      fprintf(f,"k+=%d;\n", SDF[i0][3]);

      fprintf(f,"}\n");

    }

  // this version of SDF assumes symmetry, the following will handle the upper half of the matrix.
  if(SDF[NNZD-1][0]!=0){
      //base (first NNZ)
      fprintf(f," k=%d;\n", SDF[NNZD-1][1]);

      fprintf(f," for(int j=0; j<%d; j++){\n", N/SDF[NNZD-1][3]);

      // body of N nested loops where N=depth
      for (int i1=0; i1<SDF[NNZD-1][2]-1; i1++)
	{
	  fprintf(f," for (int i%d=0; i%d<%d; i%d++){\n", i1, i1, SDF[NNZD-1][3+SDF[NNZD-1][2]+i1], i1);
	}

      // index of the next NNZ
      fprintf(f,"current_position=k");
      for (int i1=0; i1<SDF[NNZD-1][2]-1; i1++)
	{
	  fprintf(f,"+i%d*%d", i1, SDF[NNZD-1][4+i1]);
	}
      fprintf(f,";\n");

      // if not out of bound do the computation
      fprintf(f," if (current_position<%d){\n", SDF[NNZD-1][0]+N);

      fprintf(f,"C[current_position] += -%.*e * B[current_position];\n", Digs, J);

      fprintf(f," }\n");

      for (int i1=0; i1<SDF[NNZD-1][2]-1; i1++)
	{
	  fprintf(f," }\n");
	}

      // move base to next location
      fprintf(f,"k+=%d;\n", SDF[NNZD-1][3]);

      fprintf(f,"}\n");

    }

  for (int i0=0; i0<NNZD-1; i0++)
    {
      //base (first NNZ)
      fprintf(f," k=%d;\n", SDF[i0][1]);

      fprintf(f," for(int j=0; j<%d; j++){\n", N/SDF[i0][3]);

      // body of N nested loops where N=depth
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," for (int i%d=0; i%d<%d; i%d++){\n", i1, i1, SDF[i0][3+SDF[i0][2]+i1], i1);
	}

      // index of the next NNZ
      fprintf(f,"current_position=k");
      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f,"+i%d*%d", i1, SDF[i0][4+i1]);
	}
      fprintf(f,";\n");

      // if not out of bound do the computation
      fprintf(f," if (current_position<%d){\n", SDF[i0][0]+N);

      fprintf(f,"C[current_position] += -%.*e * B[current_position];\n", Digs, J);

      fprintf(f," }\n");

      for (int i1=0; i1<SDF[i0][2]-1; i1++)
	{
	  fprintf(f," }\n");
	}

      // move base to next location
      fprintf(f,"k+=%d;\n", SDF[i0][3]);

      fprintf(f,"}\n");

    }

  fprintf(f," return;\n");
  fprintf(f," }\n");

  fclose(f);

}
