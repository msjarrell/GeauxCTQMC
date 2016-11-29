#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

//#include "mpiprocess.h"

#include "init.h"
#include "det.h"


void
mpiprocess (void)
{
  printf ("mpi_wrapper: hello\n");

  int rank, size;

  MPI_Init (NULL, NULL);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  printf ("hello from %02d/%02d\n", rank, size);





  printf ("MPI finalize begin rank:%d\n", rank);
  MPI_Finalize ();

  printf ("MPI finalize finish\n");
  return;
}
