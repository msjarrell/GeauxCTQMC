#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <cuda_runtime.h>
#include <cublas_v2.h>
//#include <helper_cuda.h>

#include "det.h"


double
det (int N, double *m)
{
  cublasStatus_t status;
  cublasHandle_t handle;

  double *result_h = (double *) malloc (N * N * sizeof (double));
  int *pivot_h = (int *) malloc (N * N * sizeof (int));


  int *pivot_d;
  int *info_d;
  double *m_d;
  double *Aarray[1];


  double determinant = 1;

  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      pivot_h[i * N + j] = 0;
    }
    pivot_h[i * N + i] = 1;
  }

  //status=cublasInit();
  status = cublasCreate (&handle);
  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf (stderr, "!!!! CUBLAS initialization error\n");
    return EXIT_FAILURE;
  }


  if (cudaMalloc ((void **) &m_d, N * N * sizeof (double)) != cudaSuccess) {
    fprintf (stderr, "!!!! device memory allocation error (allocate m_d)\n");
    return EXIT_FAILURE;
  }

  Aarray[0] = m_d;

  if (cudaMalloc ((void **) &pivot_d, N * N * sizeof (int)) != cudaSuccess) {
    fprintf (stderr,
	     "!!!! device memory allocation error (allocate pivot_d)\n");
    return EXIT_FAILURE;
  }

  if (cudaMalloc ((void **) &info_d, N * N * sizeof (int)) != cudaSuccess) {
    fprintf (stderr,
	     "!!!! device memory allocation error (allocate info_d)\n");
    return EXIT_FAILURE;
  }


  printf ("Allocation Successful.\n");

  cublasSetMatrix (N, N, sizeof (double), m, N, m_d, N);

  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf (stderr, "!!!! device access error (write m_d)\n");
    return EXIT_FAILURE;
  }

  cublasSetMatrix (N, N, sizeof (double), pivot_h, N, pivot_d, N);
  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf (stderr, "!!!! device access error (write pivot_d)\n");
    return EXIT_FAILURE;
  }


  printf ("Copy to device Successful.\n");

  status = cublasDgetrfBatched (handle, N, Aarray, N, pivot_d, info_d, 1);

  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf (stderr, "!!!! kernel execution error.\n");
    return EXIT_FAILURE;
  }


  status = cublasGetMatrix (N, N, sizeof (double), Aarray[0], N, result_h, N);


  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf (stderr, "!!!! device access error (read result)\n");
    return EXIT_FAILURE;
  }


  if (cudaFree (m_d) != cudaSuccess) {
    fprintf (stderr, "!!!! memory free error (m)\n");
    return EXIT_FAILURE;
  }
  if (cudaFree (pivot_d) != cudaSuccess) {
    fprintf (stderr, "!!!! memory free error (pivot)\n");
    return EXIT_FAILURE;
  }
  if (cudaFree (info_d) != cudaSuccess) {
    fprintf (stderr, "!!!! memory free error (info)\n");
    return EXIT_FAILURE;
  }


  status = cublasDestroy (handle);



  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf (stderr, "!!!! shutdown error\n");
    return EXIT_FAILURE;
  }


  for (i = 0; i < N; i++)
    determinant *= result_h[i * N + i];

  return determinant;

}
