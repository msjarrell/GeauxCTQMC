#include "matrix.h"
#include "parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "main.h"


#define NUM 228
#define TEST_ITER 1000
#define TILE 10
//fwds 

//double random();

unsigned int seed[NUM];

void det_update_fast(Matrix *m, Matrix *temp, double *u,double *v);

double myrand(unsigned int *myseed){
  return ((double) rand_r(myseed)) / RAND_MAX;
}


void det_update_fast(Matrix *m, Matrix *temp, double *u,double *v){
  //Do the update and calculate determinat according to sherman-morrison formula
  //Assuming m->mat is an non-empty matrix, and temp->mat will be the same size
  //do temp->mat=m->mat + u X v, update inverse and determinant
  int i,j,ii,jj;
  int n=m->N;
  int sz=m->max_sz;
  double *a=m->mat;
  double *ai=m->inv;
  double *nm=temp->mat;
  double *ni=temp->inv;
  //double *aiu = malloc(sizeof(double)*n);
  //double *vtai = malloc(sizeof(double)*n);
  double aiu[800] __attribute__((aligned(64)));
  double vtai[800] __attribute__((aligned(64)));
    
  double lambda = 0 ;


  
  __assume_aligned(u,64);
  __assume_aligned(v,64);
  //__assume_aligned(a,64);
  //__assume_aligned(ai,64);
  //__assume_aligned(nm,64);
  //__assume_aligned(ni,64);  

  for(i=0;i<n;i++){
    aiu[i] = 0 ;
    vtai[i] = 0 ;
  }

  for(j=0;j<n;j++){
    double uj=u[j];
    #pragma simd 
    for(i=0;i<n;i++)
      aiu[i] += ai[i*sz+j] * uj;
      
    
    double vj=v[j];
    #pragma simd
    //#pragma vector aligned  
    for(i=0;i<n;i++)
      vtai[i] += vj * ai[j*sz+i];
  }


  for(i=0;i<n;i++){
    double ui=u[i];
#pragma simd
    //#pragma vector aligned  
    for (j=0;j<n;j++)
      nm[i*sz+j] = a[i*sz+j] + ui * v[j];
  }
  
  for(i=0;i<n;i++)
    lambda += aiu[i]*v[i];



  for(i=0;i<n;i++){
    double aiui=aiu[i];
#pragma simd 
    //#pragma vector aligned  
    for(j=0;j<n;j++)
      ni[i*sz+j]=ai[i*sz+j]-aiui*vtai[j]/(1.0+lambda);
  }

  temp->det= m->det * (1.0+lambda);

  //free(aiu);
  //free(vtai);
}


void main(){
  int task;
  const int sz=800;
  size_t msz=sizeof(double) * sz*sz;
  
  for(task=0;task<NUM;task++)
    seed[task]=time(NULL) ^ task;

  Matrix m_array[NUM];
  Matrix temp_array[NUM];


  double u_array[NUM][sz] __attribute__((aligned(64)));
  double v_array[NUM][sz] __attribute__((aligned(64)));
  /*
  double matrix_mat[NUM][msz] __attribute__((aligned(64)));
  double matrix_inv[NUM][msz] __attribute__((aligned(64)));
  double temp_mat[NUM][msz] __attribute__((aligned(64)));
  double temp_inv[NUM][msz] __attribute__((aligned(64)));
  */
  

  for(task=0;task<NUM;task++){

    Matrix *matrix=m_array+task;
    Matrix *temp=temp_array+task;
    /*
    matrix->inv=matrix_inv[task];
    matrix->mat=matrix_mat[task];
    */
    //matrix->inv=(double *)_mm_malloc(msz,64);
    matrix->inv=(double *)malloc(msz);
    //matrix->mat=(double *)_mm_malloc(msz,64);
    matrix->mat=(double *)malloc(msz);
    matrix->N=sz;
    matrix->max_sz=sz;
    matrix->det=1;
    /*
    temp->inv=temp_inv[task];
    temp->mat=temp_mat[task];
    */
    //temp->inv=(double *)_mm_malloc(msz,64);
    temp->inv=(double *)malloc(msz);
    //temp->mat=(double *)_mm_malloc(msz,64);
    temp->mat=(double *)malloc(msz);
    //m_array[task]=matrix;
    //temp_array[task]=temp;

    //u_array[task]=malloc(sizeof(double)*sz);
    //v_array[task]=malloc(sizeof(double)*sz);

  }



#pragma omp parallel for private(task)
  for(task=0;task<NUM;task++){
    printf("task %d starting\n",task);
   
    Matrix *matrix=m_array+task;
    Matrix *temp=temp_array+task;
    
    
    int i;
    for (i=0;i<sz*sz;i++){
      matrix->mat[i]=myrand(seed+task);
      matrix->inv[i]=myrand(seed+task);
    }

    int iter;
    double *u,*v;
    u=u_array[task];
    v=v_array[task];

    printf("task %d initialized\n",task);

    for(iter=0;iter<TEST_ITER;iter++){
      for(i=0;i++;i<sz){
	u[i]=myrand(seed+task);
	v[i]=myrand(seed+task);
      }

     det_update_fast(matrix,temp,u,v);
    
    }

    printf("task %d ends\n",task);
  }


  printf("cleanup...\n");
  
    for(task=0;task<NUM;task++){
      Matrix *matrix=m_array+task;
      Matrix *temp=temp_array+task;
      /*
      _mm_free(matrix->mat);
      _mm_free(matrix->inv);
      _mm_free(temp->inv);
      _mm_free(temp->mat);
      */
      free(matrix->mat);
      free(matrix->inv);
      free(temp->inv);
      free(temp->mat);


      //free(matrix);
      //free(temp);
    //free(u_array[task]);
    //free(v_array[task]);

    //free(m_array);
    //free(temp_array);
    //free(u_array);
    //free(v_array);

    }
  
    printf("success\n");
}
