#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include <math.h>

#include "matrix.h"
#include "parameters.h"

#include "det.h"
#include "utility.h"
#include "calc_green_func.h"
#include "pcg-c-basic-0.9/pcg_basic.h"


#ifdef USE_PCG
double myrand(pcg32_random_t *myseed){
  return ldexp(pcg32_random_r(myseed), -32);
}
#else
double myrand(unsigned int *myseed){
  return ((double) rand_r(myseed)) / RAND_MAX;
}
#endif


double max(double a,double b){
  return (a>b)?a:b;
}

double min(double a,double b){
  return (a<b)?a:b;
}


void print_seg(Mempool *mem){
  int i,n;
  int io;
  for(io=0;io<mem->n_orbit*2;io++){
    n=mem->n[io];
    if(n> ARRAY_SZ){
      printf("n overflow! %d\n",n);
      exit(0);
    }
    printf("%d channel:\n",io);
    for(i=0;i<n;i++)
      printf("%g\t\t%g\n",mem->is[io][i],mem->ie[io][i]);
    
  }
  printf("----------------\n");
}


void print_matrix(Matrix *m){
  int i,j;
  int n=m->N;
  printf("[");
  for(i=0;i<n;i++){
    printf("[");
    for(j=0;j<n;j++){
      int k=i*n+j;
      printf("%7g,",creal(m->mat[k]));
    }
    printf("\b],");
  }
    printf("\b]");
  printf("\ndet=%g\n",creal(m->det));
  printf("---------------------------\n");
  printf("[");
  for(i=0;i<n;i++){
    printf("[");
    for(j=0;j<n;j++){
      int k=i*n+j;
      printf("%7g,",creal(m->inv[k]));
    }
    printf("\b],");
  }
  printf("\b]\n");
  printf("---------------------------\n");
  printf("---------------------------\n");
}


void
propose_time (SEEDTYPE *seed,double *result, double left, double right, double beta)
{
  //propose a time between left and right, considering periodic boundary
  double r = myrand(seed);
  double newtime = 0.0;

  if (right <= left)
    right += beta;

  newtime = left + (right - left) * r;

  if (newtime >= beta)
    newtime -= beta;

  //printf("left=%lf,right=%lf,beta=%lf,newtime=%lf\n",left,right,beta,newtime);


  *result = newtime;
}


double overlap(double as,double ae,double bs,double be){
  double start = max(as,bs);
  double end = min(ae,be);
  return max(0.0,(end-start));
}

double overlap_total ( double *ts, double *te, int nt, double start,
		       double end, double beta){

  //printf("%g\t%g\n",start,end);
  if (start>end)
    return (overlap_total(ts,te,nt,0.0,end,beta) + overlap_total(ts,te,nt,start,beta,beta));
  else{
    int i;
    double ol=0.0;
    for (i=0;i<nt;i++){
      double as=ts[i];
      double ae=te[i];
      //  printf("%d\t%g\t%g\n",i,as,ae);
      if(ae>as)
	ol += overlap(as,ae,start,end);
      else
	ol += overlap(0.0,ae,start,end) + overlap(as,beta,start,end);
    }
    //printf("ol=%g\n",ol);
    return ol;
  } 
}


double length_total(double *ts, double *te,int nt, double beta){
    int i;
    double ol=0;
    for (i=0;i<nt;i++){
      double as=ts[i];
      double ae=te[i];
      if(ae>as)
	ol+= ae-as;
      else
	ol+= ae-as+beta;
    }
    return ol;
}

/*
void
overlap_single_single (double *lov, double is, double ie, double start,
		       double end)
{
  //calculating the overlap of two segments, most simple cases
  double overlap = 0;

  if ((start >= ie) || (is >= end)) {
    //        S----E
    //IS---IE
    //or
    //S---E
    //        IS---IE
    overlap=0;
  }
  else if (is <= start && ie >= end) {
  //    S--E
  // IS------IE
    overlap = end - start;
  }
  else if (is <= start && ie <= end) {
  //      S------E
  // IS------IE
    overlap = ie - start;
  }
  else if (start <= is && end >= ie) {
  // S------------E
  //    IS---IE
    overlap = ie - is;
  }
  else if (start <= is && end <= ie) {
  // S------E
  //    IS-----IE
    overlap = end - is;

  }

  *lov = overlap;

}





void
overlap_single (double *lov, double *ts, double *te, int nt, double start,
		double end, double beta)
{
  //overlap of single segment and all segments on the OTHER channel
  int i;
  double overlap = 0.0;

  if (end < start)
    end += beta;
  
  //printf("\nnew segment %g\t%g\n",start,end);

  

  for (i = 0; i < nt; i++) {
    
    double iov = 0;
    double is = ts[i];
    double ie = te[i];
    
    //printf("%d:%g\t%g\n",i,is,ie);

    if (ie > is) {
      overlap_single_single (&iov, is, ie, start, end);
      //printf("overlap with %dth segment %g\t%g:%g\n",i,is,ie,iov);
      overlap += iov;
      overlap_single_single (&iov, is+beta, ie+beta, start, end);
      overlap += iov;
      //printf("overlap with %dth segment %g\t%g:%g\n",i,is+beta,ie+beta,iov);
    }
    else {
      ie += beta;
      overlap_single_single (&iov, is, ie, start, end);
      //printf("overlap with %dth segment %g\t%g:%g\n",i,is,ie,iov);
      overlap += iov;
      overlap_single_single (&iov, is-beta, ie-beta, start, end);
      //printf("overlap with %dth segment %g\t%g:%g\n",i,is-beta,ie-beta,iov);

      overlap += iov;
    }
   
  }
  //  printf("%g\n",overlap);
  *lov = overlap;
}
*/

void cvt_to_MKL(MTYPE *a,MKLTYPE *m,int N){
  int i;
  
  #pragma simd  
  for(i = 0; i < N*N ; i ++ ){
    m[i].real=creal(a[i]);
    m[i].imag=cimag(a[i]);
  }
  
}


void cvt_from_MKL(MTYPE *a,MKLTYPE *m,int N){
  int i;

  #pragma simd  
  for(i=0;i<N*N;i++){
    a[i]=m[i].real+ I *m[i].imag;
  }

}



double inverse(Matrix * m)
{
  
  //  printf("Inverse\n");
  
  MTYPE *A=m->mat;
  MTYPE *B=m->inv;

  MKLTYPE *M=m->mkl;

  int N=m->N;
  int i;
  MTYPE det=1.0;

  if(N==0)
    return (0);
  
  int IPIV[ARRAY_SZ] ;
  
  cvt_to_MKL(A,M,N);
  
  LAPACKE_zgetrf( LAPACK_ROW_MAJOR, N, N, M, N, IPIV);
  
  cvt_from_MKL(B,M,N);

  for (i=0;i<N;i++)
    det*=B[i*N+i];
    
  LAPACKE_zgetri(LAPACK_ROW_MAJOR,N,M,N,IPIV);
  
  cvt_from_MKL(B,M,N);
  double error=(cabs(det)-cabs(m->det))/cabs(det);
  m->det=det;

  return error;
}


void print_list(int n, int * o_list,double *t_list){
  int i;
  printf("#=%d\n",n);
  for(i=0;i<n;i++){
    printf("%d:\t%d\t%lf\n",i,o_list[i],t_list[i]);
  }
  printf("End of list\n");
}


void verify_matrix( int i_o, Mempool *mem, Par par){
  print_seg(mem);
  double * ts=mem->is[i_o];
  double * te=mem->ie[i_o];
  int n = mem->n[i_o];
  int i,j;
  construct_matrix(ts,te,n,mem->m_tmp2,par);
  inverse(mem->m_tmp2);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      if (cabs((mem->m[i_o]->mat[i*n+j] - mem->m_tmp2->mat[i*n+j])/mem->m_tmp2->mat[i*n+j]) > 1e-3){
	printf("ERROR: matrix is not cosistent\n");
	exit(0);
      }
      if (cabs((mem->m[i_o]->inv[i*n+j] - mem->m_tmp2->inv[i*n+j])/mem->m_tmp2->inv[i*n+j]) > 1e-3){
	printf("ERROR: inverse matrix is not cosistent\n");
	exit(0);
      }
    }

}

