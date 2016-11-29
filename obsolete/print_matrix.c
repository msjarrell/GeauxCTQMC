#include <stdio.h>

#include "print_matrix.h"


void
print_matrix (int N, double *m)
{
  int i, j;
  printf ("The matrix is %d by %d:\n", N, N);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf ("%8.3f\t", m[i * N + j]);
    printf ("\n");
  }


}
