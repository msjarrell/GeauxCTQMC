#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "matrix.h"


//void construct_matrix (double **is, double **ie, int *n, Matrix * ma,MTYPE * g_tau, double beta);
void destroy_matrix (Matrix * ma);
void mv_matrix (Matrix * md, Matrix * ms);
void new_matrix (Matrix *m, int n);
void copy_matrix (Matrix * md, Matrix * ms);

void new_matrix(Matrix *m,int n){
  m->max_sz=n;
  m->N=n;
  m->det=1.0;
  m->mat=malloc(sizeof(MTYPE) *n *n);
  m->inv=malloc(sizeof(MTYPE) *n *n);
}


void
destroy_matrix (Matrix * m)
{
  if (m->mat) free (m->mat);
  m->mat = NULL;
  if(m->inv) free (m->inv);
  m->inv = NULL;
  m->max_sz = 0;
  m->N = 0;
  m->det = 0;
}

void
mv_matrix (Matrix * md, Matrix * ms)
{
  if (md->mat) free (md->mat);
  if (md->inv) free (md->inv);
  md->mat = ms->mat;
  md->inv = ms->inv;
  md->max_sz = ms->max_sz;
  md->N = ms->N;
  md->det = ms->det;
}


void
copy_matrix (Matrix * md, Matrix * ms)
{
  int i,j;
  int sz=ms->max_sz;
  size_t size=sz*sz*sizeof(MTYPE);
  md->mat=malloc(size);
  md->inv=malloc(size);

  memcpy(md->mat,ms->mat,size);
  memcpy(md->inv,ms->inv,size);
  /*
  for(i=0;i<sz;i++)
    for(j=0;j<sz;j++){
      md->mat[i*sz+j] = ms->mat[i*sz+j];
      md->inv[i*sz+j] = ms->inv[i*sz+j];
    }
  */
  md->max_sz = ms->max_sz;
  md->N = ms->N;
  md->det = ms->det;
}

/*
void
construct_matrix (double **is, double **ie, int *n, Matrix * m,
		  MTYPE * g_tau, double beta)
{

  int i, j;

  int sz = (INIT_SZ > *n) ? INIT_SZ : (*n);
  m->max_sz = sz;


  m->mat = malloc (sizeof (MTYPE) * sz * sz);
  m->inv = malloc (sizeof (MTYPE) * sz * sz);

  m->N = *n;
  for (i = 0; i < sz * sz; i++) {
    m->mat[i] = 0 + 0 * I;
    m->inv[i] = 0 + 0 * I;
  }
  //put in the matrix elements.
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      double start = is[0][i];
      double end = ie[0][j];
      double tau = end - start;
      MTYPE result = 0;
      get_g_tau (g_tau, tau, beta, &result);
      m->mat[i * sz + j] = result;
    }
  }
  //calculate determinant
  det (m->max_sz, m->N, m->mat, &(m->det));

}
*/
