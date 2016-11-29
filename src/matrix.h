#include <complex.h>    /* Standart Library of Complex Numbers */
#include <mkl.h>
#include "parameters.h"

#ifdef USE_PCG
#include "pcg-c-basic-0.9/pcg_basic.h"
#endif

//Default matrix size
#define INIT_SZ 1

//#define rand_r(a) 100


#ifdef USE_PCG
typedef pcg32_random_t SEEDTYPE;
#else
typedef unsigned int SEEDTYPE;
#endif

typedef double complex MTYPE ;
typedef MKL_Complex16  MKLTYPE;
//typedef double MTYPE ;

//typedef double SDL_MATRIX_TYPE ;


typedef struct
{
  double V[HM_SZ];//diagonal matrix with eigenvalues of the hamiltonian
  double U[HM_SZ*HM_SZ];//eigenvectors of H 
  double OP[4][HM_SZ*HM_SZ];//each of OP matrix corresponds to a U'FU, F is a fermion operator matrix
  
  int proc_id;
  int job_id;
  double beta;
  double u;
  double ef;
  double v2;
  double D;
  MTYPE g_tau[N_TAU];
  int iter_measure;
  int check_err;
  int n_check;
  int iter_warm;
  double t0;
  double t1;
  int task;
  int seed;
#ifdef SDF
  double vecd[HM_SZ];
#endif
} Par;

typedef struct
{
  MTYPE *mat;//The matrix itself
  int max_sz;//Physical size for storage
  int N;//logical size, N<=max_sz
  MTYPE *inv;//The inverse of the matrix
  MTYPE det;//The determinant of the matrix
  MKLTYPE *mkl;//storage for MKL-compatible data
} Matrix;

typedef struct
{
  //SDL_MATRIX_TYPE *sdl_m

#ifdef USE_MATRIX
  double trace;
#endif

  int n_orbit;

  Matrix *m[N_ORBIT*2];
  Matrix *m_tmp1;
  Matrix *m_tmp2;

  //number of operator pairs on each orbital
  int n[N_CHANNEL];
  int isFull[N_CHANNEL];
  double n_av[N_CHANNEL];
  double err_max;
  double nn_av[N_CHANNEL];

  //total number of operators on all orbitals, n_sum=2 * sum(n[])
  int n_sum;

  //list of operators used in diagonal models
  double __attribute__((aligned(64))) is[N_CHANNEL][ARRAY_SZ];
  double __attribute__((aligned(64))) ie[N_CHANNEL][ARRAY_SZ];



#ifdef USE_MATRIX
  //list of operators used in non-diagonal models
  int *op_list;
  double *time_list;
  int *tmp_op_list;
  double *tmp_time_list;
#endif
  //temp matrix and vector for local term evaluation used in non-diagonal models
  //  double *m_loc;
  //double *evt;
  
  //matrices for using MKL dgemm
  //  double *MKL_temp;//, *MKL_B, *MKL_C;

  //measured quantities

#ifdef MEASURE_GIO
  double complex green_iomega[N_CHANNEL][2*N_OMEGA];
#endif

  double complex green_tau[N_CHANNEL][N_TAU];
  unsigned int green_tau_n[N_CHANNEL][N_TAU];

#ifdef MEASURE_SCPT
  double chargeScpt[N_ORBIT][N_TAU];
  double spinScpt[N_ORBIT][N_TAU];
#endif

#ifdef MEASURE_LEG
  double complex green_legendre[N_CHANNEL][N_LEG] ;
#endif
  //utility variables

  int id;
  int iter;

  SEEDTYPE seed;

  //counters for statistics
  int add_prop;
  int add_accpt;
  int remove_prop;
  int remove_accpt;
  int shift_prop;
  int shift_accpt;
  int add_anti_prop;
  int add_anti_accpt;
  int rem_anti_prop;
  int rem_anti_accpt;

  int add_11;
  int add_10;
  int add_01;
  int add_00;

  int rem_11;
  int rem_10;
  int rem_01;
  int rem_00;

} Mempool;
