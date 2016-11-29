#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#include "matrix.h"



void copy_matrix (Matrix * m, Matrix * n);
void add_row (double **ts, double **te, int nt, int newid, double news,
	      double newe, Matrix * m);

void
copy_matrix (Matrix * m, Matrix * n)
{

  n->mat = malloc (sizeof (MTYPE) * INIT_SZ * INIT_SZ);
  n->inv = malloc (sizeof (MTYPE) * INIT_SZ * INIT_SZ);

  n->max_sz = m->max_sz;
  n->N = m->N;
  n->det = m->det;

  memcpy (n->mat, m->mat, sizeof (MTYPE) * INIT_SZ * INIT_SZ);
}


void
expand_matrix (Matrix * m)
{

  int new_size = m->max_sz + INIT_SZ;
  MTYPE *new_space =
    malloc (sizeof (MTYPE) * new_size * new_size);
  MTYPE *new_space_inv =
    malloc (sizeof (MTYPE) * new_size * new_size);

  int i, j;
  for (i = 0; i < m->N; i++)
    for (j = 0; j < m->N; j++) {
      new_space[i * new_size + j] = m->mat[i * m->max_sz + j];
      new_space_inv[i * new_size + j] = m->inv[i * m->max_sz + j];
    }
  m->max_sz = new_size;
  free (m->mat);
  free (m->inv);
  m->mat = new_space;
  m->inv = new_space_inv;
}



void
add_row (double **ts, double **te, int nt, int newid, double news,
	 double newe, Matrix * m)
{
  if (nt == m->max_sz)
    expand_matrix (m);

  int n = m->N;
  int SZ = m->max_sz;
  int i, j;

  for (i = 0; i < n; i++)
    for (j = n - 1; j >= newid; j++)
      m->mat[i * SZ + j + 1] = m->mat[i * SZ + j];

  for (j = 0; j < n + 1; j++)
    for (i = n - 1; i >= newid; i++)
      m->mat[(i + 1) * SZ + j] = m->mat[i * SZ + j];


  for (i = 0; i < newid; i++) {
    MTYPE result = 0;
    double tau = fabs (ts[0][i] - newe);
    get_g_tau (g_tau, tau, beta, &result);
    m->mat[i * SZ + newid] = result;
  }

  for (i = newid + 1; i < n + 1; i++) {
    MTYPE result = 0;
    double tau = fabs (ts[0][i - 1] - newe);
    get_g_tau (g_tau, tau, beta, &result);
    m->mat[i * SZ + newid] = result;

  }


  for (j = 0; j < new_id; j++) {
    MTYPE result = 0;
    double tau = fabs (te[0][j] - news);
    get_g_tau (g_tau, tau, beta, &result);
    m->mat[j + newid * SZ] = result;
  }

  for (j = new_id + 1; j < n + 1; j++) {
    MTYPE result = 0;
    double tau = fabs (te[0][j - 1] - news);
    get_g_tau (g_tau, tau, beta, &result);
    m->mat[j + newid * SZ] = result;

  }
  MTYPE result = 0;
  double tau = fabs (news - newe);
  get_g_tau (g_tau, tau, beta, &result);
  m->mat[newid * SZ + new_id] = result;

  m->N = n + 1;
  det (m->max_sz, m->N, m->mat, &(m->det));

}


void
del_row (double **ts, double **te, int nt, int remid, Matrix * m)
{

  int n = m->N;
  int SZ = m->max_sz;
  int i, j;

  for (i = 0; i < n; i++)
    for (j = remid; j < n; j++)
      m->mat[i * SZ + j] = m->mat[i * SZ + j + 1];

  for (j = 0; j < n; j++)
    for (i = remid; i < n; i++)
      m->mat[i * SZ + j] = m->mat[(i + 1) * SZ + j];

  m->N = n - 1;
  det (m->max_sz, m->N, m->mat, &(m->det));

}

void
mod_row (double **ts, double **te, int nt, int modid, double news,
	 double newe, Matrix * m)
{

  int n = m->N;
  int SZ = m->max_sz;
  int i, j;

  for (i = 0; i < n; i++) {
    MTYPE result = 0;
    double tau = fabs (ts[0][i] - newe);
    get_g_tau (g_tau, tau, beta, &result);
    m->mat[i * SZ + modid] = result;
  }



  for (i = 0; i < n; i++) {
    MTYPE result = 0;
    double tau = fabs (news - te[0][i]);
    get_g_tau (g_tau, tau, beta, &result);
    m->mat[i + SZ * modid] = result;
  }

  MTYPE result = 0;
  double tau = fabs (news - newe);
  get_g_tau (g_tau, tau, beta, &result);
  m->mat[modid + SZ * modid] = result;


  det (m->max_sz, m->N, m->mat, &(m->det));

}
