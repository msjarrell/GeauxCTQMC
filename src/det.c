#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>

#include "matrix.h"
#include "parameters.h"


#include "det.h"
#include "calc_green_func.h"

// For these methods; storage[0][0] is the determinant before the change is made, 
//                      storage[0][1] is the determinant after the change is made


//ts = temp interval start times
//te = temp interval end times
//nt = number of temp intervals
//up = 1 if the intervals are from the UP spins, 0 if from DOWN

//newid = id of value being added (3 = before previous interval with id=3)
//news = new start time
//newe = new end time

void det_update_fast(Matrix *m, Matrix *temp, MTYPE *u,MTYPE *v){
  //Do the update and calculate determinat according to sherman-morrison formula
  //Assuming m->mat is an non-empty matrix, and temp->mat will be the same size
  //do temp->mat=m->mat + u X v, update inverse and determinant
  int i,j;
  int n=m->N;
  //int sz=m->max_sz;
  int sz=n;
  temp->N=n;
  //size_t array_sz=sizeof(MTYPE) *n;
  MTYPE *a=m->mat;
  MTYPE *ai=m->inv;
  MTYPE *nm=temp->mat;
  MTYPE *ni=temp->inv;
  //MTYPE *aiu = (MTYPE *)_mm_malloc(array_sz,64);
  //MTYPE *vtai = (MTYPE *)_mm_malloc(array_sz,64);

  assert(n<ARRAY_SZ);
  MTYPE aiu[ARRAY_SZ] __attribute__((aligned(64)));
  MTYPE vtai[ARRAY_SZ] __attribute__((aligned(64)));
  
  
  __assume_aligned(u,64);
  __assume_aligned(v,64);
  __assume_aligned(a,64);
  __assume_aligned(ai,64);
  __assume_aligned(nm,64);
  __assume_aligned(ni,64);;
  //__assume_aligned(aiu,64);
  //__assume_aligned(vtai,64);

  MTYPE lambda = 0 ;
  
  for(i=0;i<n;i++){
    aiu[i] = 0 ;
    vtai[i] = 0 ;
  }
  
  /*

#pragma simd   
  for (i=0;i<n*n;i++){
    nm[i]=a[i];
  }

#pragma simd 
  for (i=0;i<n*n;i++){
    ni[i]=ai[i];
  }


  CBLAS_ORDER     order;
  CBLAS_TRANSPOSE trans1,trans2;

  order = CblasRowMajor;
  trans1 = CblasNoTrans;
  trans2 = CblasTrans;

  cblas_dgemv (order,trans1,sz, sz, 1.0, ai, sz, u, 1, 0, aiu, 1);
  cblas_dger(order,sz,sz,1.0,u,1,v,1,nm,sz);
  cblas_dgemv (order,trans2,sz, sz, 1.0, ai, sz, v, 1, 0, vtai, 1);
  lambda = cblas_ddot(sz,aiu,1,v,1);
  cblas_dger(order,sz,sz,-1.0/(1.0+lambda),aiu,1,vtai,1,ni,sz);
  */
    
  for(j=0;j<n;j++){
    MTYPE uj=u[j];
#pragma simd 
    for(i=0;i<n;i++)
      aiu[i] += ai[i*sz+j] * uj;
      
    
    MTYPE vj=v[j];
#pragma simd
    for(i=0;i<n;i++)
      vtai[i] += vj * ai[j*sz+i];

  }

  for(i=0;i<n;i++){
    MTYPE ui=u[i];
#pragma simd
    for(j=0;j<n;j++){
      nm[i*sz+j] = a[i*sz+j] + ui * v[j];
    }
  }


  for(i=0;i<n;i++)
    lambda += aiu[i]*v[i];

  for(i=0;i<n;i++){
    MTYPE aiui=aiu[i];
    #pragma simd 
    for(j=0;j<n;j++)
      ni[i*sz+j]=ai[i*sz+j]-aiui*vtai[j]/(1.0+lambda);


  }
  
  temp->det= m->det * (1.0+lambda);

  //free(aiu);
  //free(vtai);
  //_mm_free(aiu);
  //_mm_free(vtai);


}




void
det_add_fast (MTYPE *storage, double *ts, double *te, int nt, int newid,
	      double news, double newe, Matrix * m, Matrix * temp1, Matrix *temp2,
	 Par par)
{
 
  if (nt==0){
    //Old matrix is empty, construct new 1x1 matrix
    storage[0] = 1.0;
    temp1->max_sz=ARRAY_SZ;
    temp1->N=1;
    MTYPE result;
    double tau=newe-news;
    get_g_tau (par.g_tau, tau, par.beta, &result);
    temp1->mat[0] = result;
    temp1->inv[0] = 1.0/result;
    temp1->det = result;
    storage[1] = result;
    
  }
  else{
    //old matrix is not empty
    storage[0] = m->det;

    //construct space for A^(N+1) and inverse
    int i,j;
    int nsz=nt+1;
    //int sz=ARRAY_SZ;
    int sz=nt;
    //target matrix
    MTYPE *nm=temp1->mat;
    MTYPE *ni=temp1->inv;
    
    //temp matrix
    MTYPE *a=temp2->mat;
    MTYPE *ai=temp2->inv;

    //original matrix
    MTYPE *mm=m->mat;
    MTYPE *mi=m->inv;

    __assume_aligned(mm,64);
    __assume_aligned(mi,64);
    __assume_aligned(a,64);
    __assume_aligned(ai,64);
    __assume_aligned(nm,64);
    __assume_aligned(ni,64);


    /************************
      | A B |   |A 0 B|
      | C D |   |0 1 0|
                |C 0 D|


                |A B 0|
                |C D 0|
                |0 0 1|


    *************************/



    //part A
      for (i=0;i<newid;i++){
#pragma simd
	for (j=0;j<newid;j++)
	  nm[i*nsz+j]=mm[i*sz+j];
	
#pragma simd
	for (j=0;j<newid;j++)
	  ni[i*nsz+j]=mi[i*sz+j];
	
	nm[i*nsz+newid]=0;//right to A
	ni[i*nsz+newid]=0;
	
      //part B
#pragma simd
	for (j=newid;j<nt;j++)
	  nm[i*nsz+j+1]=mm[i*sz+j];
	
#pragma simd
	for (j=newid;j<nt;j++)
	  ni[i*nsz+j+1]=mi[i*sz+j];      
      }    
    
      
      for (i=newid;i<nt;i++){
    //part C
#pragma simd
	for (j=0;j<newid;j++)
	  nm[(i+1)*nsz+j]=mm[i*sz+j];
#pragma simd
	for (j=0;j<newid;j++)
	  ni[(i+1)*nsz+j]=mi[i*sz+j];
	
	nm[(i+1)*nsz+newid]=0;//right to C
	ni[(i+1)*nsz+newid]=0;


//part D
#pragma simd
	for (j=newid;j<nt;j++)
	  nm[(i+1)*nsz+j+1]=mm[i*sz+j];
	
#pragma simd
	for (j=newid;j<nt;j++)
	  ni[(i+1)*nsz+j+1]=mi[i*sz+j];      
      }

#pragma simd
      for (j=0;j<nsz;j++)
	nm[newid*nsz+j]=0;//Below A and B
	

#pragma simd
      for (j=0;j<nsz;j++)
	ni[newid*nsz+j]=0;//Below A and B
	

	
      nm[newid*nsz+newid]=1;
      ni[newid*nsz+newid]=1;
    



    //temp->max_sz=nsz;
    //temp->N=nsz;    
    //temp->mat=nm;
    //temp->inv=ni;
    //temp->det=m->det;

    //    Matrix *A_expand=malloc(sizeof(Matrix));

    //copy_matrix (A_expand, temp);

    //add column
    //size_t array_sz=sizeof(MTYPE) *nsz;

    //assert(nsz<ARRAY_SZ);
    MTYPE u[ARRAY_SZ] __attribute__((aligned(64)));
    MTYPE v[ARRAY_SZ] __attribute__((aligned(64)));

    //MTYPE *u=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *v=(MTYPE *)_mm_malloc(array_sz,64 );

    //MTYPE *v=malloc(sizeof(MTYPE) *nsz );
    double tau;
    MTYPE result;

    
    for (i=0;i<newid;i++){
      tau=newe-ts[i];
      get_g_tau (par.g_tau, tau, par.beta, &result);

      u[i]=result;
      v[i]=0;
    }
    for (i=newid;i<nt;i++){
      tau=newe-ts[i];
      get_g_tau (par.g_tau, tau, par.beta, &result);
      u[i+1]=result;
      v[i+1]=0;
    }
    
    tau=newe-news;
    get_g_tau (par.g_tau, tau, par.beta, &result);
    u[newid]=result-1;
    v[newid]=1;

    temp1->det=m->det;
    temp1->N=nsz;
    temp1->max_sz=ARRAY_SZ;

    det_update_fast(temp1, temp2 , u, v);//add the column, store the result in A_expand
 
    //add the row
    for (i=0;i<newid;i++){
      tau=te[i]-news;
      get_g_tau (par.g_tau, tau, par.beta, &result);
      v[i]=result;
      u[i]=0;
    }

    for (i=newid;i<nt;i++){
      tau=te[i]-news;
      get_g_tau (par.g_tau, tau, par.beta, &result);
      v[i+1]=result;
      u[i+1]=0;
    }
    
    v[newid]=0;
    u[newid]=1;

    temp2->N=nsz;
    temp2->max_sz=ARRAY_SZ;
 
    det_update_fast(temp2, temp1, u, v);//add the column, store the result in A_expand
    
    //    destroy_matrix(A_expand);
    //free(A_expand);
    //free(u);
    //free(v);

    //_mm_free(u);
    //_mm_free(v);

    storage[1]=temp1->det;
    
  }
}


void
det_add_anti_fast (MTYPE *storage, double *ts, double *te, int nt, int newid,
	      double segStart, double segEnd, Matrix * m, Matrix * temp1, Matrix *temp2,
	 Par par)
{
 
  if (nt==0){
    //Old matrix is empty, construct new 1x1 matrix
    storage[0] = 1.0;
    temp1->max_sz=ARRAY_SZ;
    temp1->N=1;
    MTYPE result;
    double news = segEnd;
    double newe = segStart;
    double tau=newe-news;
    get_g_tau (par.g_tau, tau, par.beta, &result);
    temp1->mat[0] = result;
    temp1->inv[0] = 1.0/result;
    temp1->det = result;
    storage[1] = result;
    
  }
  else{
    //old matrix is not empty
    storage[0] = m->det;

    //construct space for A^(N+1) and inverse
    int i,j;
    int nsz=nt+1;
    //int sz=ARRAY_SZ;
    int sz=nt;
    //target matrix
    MTYPE *nm=temp1->mat;
    MTYPE *ni=temp1->inv;
    
    //temp matrix
    MTYPE *a=temp2->mat;
    MTYPE *ai=temp2->inv;

    //original matrix
    MTYPE *mm=m->mat;
    MTYPE *mi=m->inv;

    __assume_aligned(mm,64);
    __assume_aligned(mi,64);
    __assume_aligned(a,64);
    __assume_aligned(ai,64);
    __assume_aligned(nm,64);
    __assume_aligned(ni,64);


    /************************
      | A B |   |A 0 B|
      | C D |   |0 1 0|
                |C 0 D|


                |A B 0|
                |C D 0|
                |0 0 1|


    *************************/

    int rowid = newid+1; 
    int colid = newid; 

    double news = segEnd;
    double newe = segStart;
    

    //part A
      for (i=0;i<rowid;i++){
#pragma simd
	for (j=0;j<colid;j++)
	  nm[i*nsz+j]=mm[i*sz+j];
	
#pragma simd
	for (j=0;j<colid;j++)
	  ni[j*nsz+i]=mi[j*sz+i];
	
	nm[i*nsz+colid]=0;//right to A
	ni[colid*nsz+i]=0;
	
      //part B
#pragma simd
	for (j=colid;j<nt;j++)
	  nm[i*nsz+j+1]=mm[i*sz+j];
	
#pragma simd
	for (j=colid;j<nt;j++)
	  ni[(j+1)*nsz+i]=mi[j*sz+i];      
      }    
    
      
      for (i=rowid;i<nt;i++){
    //part C
#pragma simd
	for (j=0;j<colid;j++)
	  nm[(i+1)*nsz+j]=mm[i*sz+j];
#pragma simd
	for (j=0;j<colid;j++)
	  ni[(j)*nsz+i+1]=mi[j*sz+i];
	
	nm[(i+1)*nsz+colid]=0;//right to C
	ni[(colid)*nsz+(i+1)]=0;


//part D
#pragma simd
	for (j=colid;j<nt;j++)
	  nm[(i+1)*nsz+j+1]=mm[i*sz+j];
	
#pragma simd
	for (j=colid;j<nt;j++)
	  ni[(j+1)*nsz+i+1]=mi[j*sz+i];      
      }

#pragma simd
      for (j=0;j<nsz;j++)
	nm[rowid*nsz+j]=0;//Below A and B
	

#pragma simd
      for (j=0;j<nsz;j++)
	ni[j*nsz+rowid]=0;//Below A and B
	

	
      nm[rowid*nsz+colid]=1;
      ni[colid*nsz+rowid]=1;
    



    //temp->max_sz=nsz;
    //temp->N=nsz;    
    //temp->mat=nm;
    //temp->inv=ni;
    //temp->det=m->det;

    //    Matrix *A_expand=malloc(sizeof(Matrix));

    //copy_matrix (A_expand, temp);

    //add column
    //size_t array_sz=sizeof(MTYPE) *nsz;

    //assert(nsz<ARRAY_SZ);
    MTYPE u[ARRAY_SZ] __attribute__((aligned(64)));
    MTYPE v[ARRAY_SZ] __attribute__((aligned(64)));

    //MTYPE *u=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *v=(MTYPE *)_mm_malloc(array_sz,64 );

    //MTYPE *v=malloc(sizeof(MTYPE) *nsz );
    double tau;
    MTYPE result;

    //add the column
    for (i=0;i<rowid;i++){
      tau=newe-ts[i];
      get_g_tau (par.g_tau, tau, par.beta, &result);

      u[i]=result;
      v[i]=0;
    }
    for (i=rowid;i<nt;i++){
      tau=newe-ts[i];
      get_g_tau (par.g_tau, tau, par.beta, &result);
      u[i+1]=result;
      v[i+1]=0;
    }
    
    tau=newe-news;
    get_g_tau (par.g_tau, tau, par.beta, &result);
    u[rowid]=result-1;
    v[rowid]=0;
    v[colid]=1;

    temp1->det=m->det;
    temp1->N=nsz;
    temp1->max_sz=ARRAY_SZ;

    det_update_fast(temp1, temp2 , u, v);//add the column, store the result in A_expand
 
    //add the row
    for (i=0;i<colid;i++){
      tau=te[i]-news;
      get_g_tau (par.g_tau, tau, par.beta, &result);
      v[i]=result;
      u[i]=0;
    }

    for (i=colid;i<nt;i++){
      tau=te[i]-news;
      get_g_tau (par.g_tau, tau, par.beta, &result);
      v[i+1]=result;
      u[i+1]=0;
    }
    
    v[colid]=0;
    u[colid]=0;
    u[rowid]=1;

    temp2->N=nsz;
    temp2->max_sz=ARRAY_SZ;
 
    det_update_fast(temp2, temp1, u, v);//add the column, store the result in A_expand
    
    //    destroy_matrix(A_expand);
    //free(A_expand);
    //free(u);
    //free(v);

    //_mm_free(u);
    //_mm_free(v);

    storage[1]=temp1->det;
    
  }
}



void det_remove_fast (MTYPE *storage, double *ts, double *te, int nt, int remid,
		      Matrix * m, Matrix * temp1,Matrix * temp2, Par par){
  storage[0] = m->det;
  //construct space for A^(N+1) and inverse
  int i,j;
  int nsz=nt-1;
  //int sz=ARRAY_SZ;
  int sz=nt;
  
  //  MTYPE * g_tau=par.g_tau;
  //  double beta=par.beta;

  if(nsz>0){  

    //remove column
    //size_t array_sz=sizeof(MTYPE) *sz;
    
    //assert(nsz<ARRAY_SZ);
    MTYPE u[ARRAY_SZ] __attribute__((aligned(64)));
    MTYPE v[ARRAY_SZ] __attribute__((aligned(64)));
    
    //MTYPE *u=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *v=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *u=malloc(sizeof(MTYPE) *sz);
    //MTYPE *v=malloc(sizeof(MTYPE) *sz);
    double tau;
    //MTYPE result;
    
    MTYPE *nm=temp1->mat;
    MTYPE *ni=temp1->inv;
  
    //temp matrix
    MTYPE *am=temp2->mat;
    MTYPE *ai=temp2->inv;

    //original matrix
    MTYPE *mm=m->mat;
    MTYPE *mi=m->inv;

    __assume_aligned(mm,64);
    __assume_aligned(mi,64);
    __assume_aligned(am,64);
    __assume_aligned(ai,64);
    __assume_aligned(nm,64);
    __assume_aligned(ni,64);
    

    //remove the id'th column  
#pragma simd
    for (i=0;i<nt;i++){
      u[i]=-mm[i*sz+remid];;
      v[i]=0;
    }
 
    u[remid]=1-mm[remid*sz+remid];
    v[remid]=1;
    
    det_update_fast(m, temp1, u, v);


    temp1->N=sz;
    temp1->max_sz=ARRAY_SZ;
    
    
    //remove the row
#pragma simd
    for (i=0;i<nt;i++){
      v[i]=-nm[remid*sz+i];
      u[i]=0;
    }
    
    v[remid]=0;
    u[remid]=1;
 
    det_update_fast(temp1, temp2, u, v);
 
    temp1->det=temp2->det;

    temp1->N=nsz;
  
    for(i=0;i<remid;i++){
#pragma simd
      for(j=0;j<remid;j++){
	nm[i*nsz+j]=am[i*sz+j];
	ni[i*nsz+j]=ai[i*sz+j];
      }
#pragma simd
      for(j=remid;j<nsz;j++){
	nm[i*nsz+j]=am[i*sz+j+1];
	ni[i*nsz+j]=ai[i*sz+j+1];
      }
    }
    
    for(i=remid;i<nsz;i++){
#pragma simd
      for(j=0;j<remid;j++){
	nm[i*nsz+j]=am[(i+1)*sz+j];
	ni[i*nsz+j]=ai[(i+1)*sz+j];
      }
#pragma simd
      for(j=remid;j<nsz;j++){
	nm[i*nsz+j]=am[(i+1)*sz+j+1];
	ni[i*nsz+j]=ai[(i+1)*sz+j+1];
      }
    }

    storage[1]=temp1->det;
  }    
  
  else{
    //new matrix is empty

    temp1->max_sz=ARRAY_SZ;
    temp1->N=0;
    temp1->det = 1.0;
    storage[1]=1.0;
  }

}



void
det_shift_fast (MTYPE *storage, double *ts, double *te, int nt, int modid,
		 double news, double newe, Matrix * m, Matrix * temp1, Matrix *temp2,
	    Par par)
{
  storage[0] = m->det;
  storage[1] = 1.0;

  int sz=nt;
  
  MTYPE u[ARRAY_SZ] __attribute__((aligned(64)));
  MTYPE v[ARRAY_SZ] __attribute__((aligned(64)));

  int i;

  double beta=par.beta;

  MTYPE *mm=m->mat;
  MTYPE *mi=m->inv;
  
  MTYPE *nm=temp1->mat;
  MTYPE *ni=temp1->inv;

  //temp matrix
  MTYPE *am=temp2->mat;
  MTYPE *ai=temp2->inv;
  
  __assume_aligned(mm,64);
  __assume_aligned(mi,64);
  __assume_aligned(nm,64);
  __assume_aligned(ni,64);

  __assume_aligned(am,64);
  __assume_aligned(ai,64);


  //change column
#pragma simd
  for (i=0;i<nt;i++){
    double tau=newe-ts[i];
    MTYPE result;
    get_g_tau (par.g_tau, tau, beta, &result);
    u[i]=result - mm[i*sz+modid];
    v[i]=0;
    }
 
  double tau=newe-news;
  MTYPE result;
  get_g_tau (par.g_tau, tau, beta, &result);
  u[modid]=result-mm[modid*sz+modid];
  v[modid]=1;

  det_update_fast(m, temp2, u, v);

  //printf("det=%g\n",m->det);
  //printf("det=%g\n",temp2->det);

  //change the row
#pragma simd
    for (i=0;i<nt;i++){
      double tau=te[i]-news;
      MTYPE result;
      get_g_tau (par.g_tau, tau, beta, &result);
      v[i]=result-am[modid*sz+i];
      u[i]=0.0;
    }
    
    v[modid]=0.0;
    u[modid]=1;
 
    det_update_fast(temp2, temp1, u, v);
    // printf("det=%g\n",temp1->det);

  storage[1] = temp1->det;
}




void
back_det_add_anti_fast (MTYPE *storage, double *ts, double *te, int nt, int newid,
	      double news, double newe, Matrix * m, Matrix * temp1, Matrix *temp2,
	 Par par)
{

  
    //old matrix is not empty
  storage[0] = m->det;

  double beta=par.beta;

  int i,j;
  int nsz=nt+1;
  //int sz=ARRAY_SZ;
  int sz=nt;
  //target matrix
  MTYPE *nm=temp1->mat;
  MTYPE *ni=temp1->inv;
  
  //temp matrix
  MTYPE *a=temp2->mat;
  MTYPE *ai=temp2->inv;
  
  //original matrix
  MTYPE *mm=m->mat;
  MTYPE *mi=m->inv;
  
  MTYPE u[ARRAY_SZ] __attribute__((aligned(64)));
  MTYPE v[ARRAY_SZ] __attribute__((aligned(64)));


  __assume_aligned(mm,64);
  __assume_aligned(mi,64);
  __assume_aligned(a,64);
  __assume_aligned(ai,64);
  __assume_aligned(nm,64);
  __assume_aligned(ni,64);

  
  //for(i=0;i<nt;i++)
  //  printf("%g\t%g\n",ts[i],te[i]);

  
  //printf("-------------\n%g\t%g\n",news,newe);
  //  printf("id=%d\n",newid);
  //  double t_end=te[newid];
  double tau;
  MTYPE result;

  //printf("Step1\n");

  //change the row
#pragma simd
  for (i=0;i<nt;i++){
    double tau=te[i]-newe;
    MTYPE result;
    get_g_tau (par.g_tau, tau, beta, &result);
    v[i]=result-mm[newid*sz+i];
    u[i]=0.0;
  }


  u[newid]=1;

  det_update_fast(m, temp2, u, v);
  //print_matrix(m);
  //  print_matrix(temp2);
  //  printf("det=%g\n",m->det);
  //printf("det=%g\n",temp2->det);





//      | A B |   |A 0 B|
//      | C D |   |0 1 0|
//                |C 0 D|


//                |A B 0|
//                |C D 0|
//                |0 0 1|




  //printf("Step2\n");

    //part A
  for (i=0;i<newid;i++){
#pragma simd
    for (j=0;j<newid;j++)
      nm[i*nsz+j]=a[i*sz+j];
    
#pragma simd
    for (j=0;j<newid;j++)
      ni[i*nsz+j]=ai[i*sz+j];
    
    nm[i*nsz+newid]=0;//right to A
    ni[i*nsz+newid]=0;
	
      //part B
#pragma simd
    for (j=newid;j<nt;j++)
      nm[i*nsz+j+1]=a[i*sz+j];
    
#pragma simd
    for (j=newid;j<nt;j++)
      ni[i*nsz+j+1]=ai[i*sz+j];      
  }    
    
      
  for (i=newid;i<nt;i++){
    //part C
#pragma simd
    for (j=0;j<newid;j++)
      nm[(i+1)*nsz+j]=a[i*sz+j];
#pragma simd
    for (j=0;j<newid;j++)
      ni[(i+1)*nsz+j]=ai[i*sz+j];
    
    nm[(i+1)*nsz+newid]=0;//right to C
    ni[(i+1)*nsz+newid]=0;


//part D
#pragma simd
    for (j=newid;j<nt;j++)
      nm[(i+1)*nsz+j+1]=a[i*sz+j];
    
#pragma simd
    for (j=newid;j<nt;j++)
      ni[(i+1)*nsz+j+1]=ai[i*sz+j];      
  }

#pragma simd
  for (j=0;j<nsz;j++)
    nm[newid*nsz+j]=0;//Below A and B
	

#pragma simd
  for (j=0;j<nsz;j++)
    ni[newid*nsz+j]=0;//Below A and B
	

	
  nm[newid*nsz+newid]=1;
  ni[newid*nsz+newid]=1;


  //printf("Step3\n");    

  for (i=0;i<newid;i++){
    tau=news-ts[i];
    get_g_tau (par.g_tau, tau, par.beta, &result);
    u[i]=result;
    v[i]=0;
  }
  for (i=newid;i<nt;i++){
    tau=news-ts[i];
    get_g_tau (par.g_tau, tau, par.beta, &result);
    u[i+1]=result;
    v[i+1]=0;
  }
    
  tau=news-newe;
  get_g_tau (par.g_tau, tau, par.beta, &result);
  u[newid+1]=result;
  

  tau=news-ts[newid];
  get_g_tau (par.g_tau, tau, par.beta, &result);
  u[newid]=result-1;
  v[newid]=1;

  temp1->det=temp2->det;
  temp1->N=nsz;
  temp1->max_sz=ARRAY_SZ;

  //  print_matrix(temp1);
  
  det_update_fast(temp1, temp2 , u, v);//add the column, store the result in A_expand
  //printf("det=%g\n",temp2->det);

  //  print_matrix(temp2);

  //printf("Step4\n"); 
    //add the row
  for (i=0;i<newid;i++){
    tau=te[i]-ts[newid];
    get_g_tau (par.g_tau, tau, par.beta, &result);
    v[i]=result;
    u[i]=0;
  }

  for (i=newid;i<nt;i++){
    tau=te[i]-ts[newid];
    get_g_tau (par.g_tau, tau, par.beta, &result);
    v[i+1]=result;
    u[i+1]=0;
  }
    
  v[newid]=0;
  u[newid]=1;
  
  
  det_update_fast(temp2, temp1, u, v);//add the column, store the result in A_expand
  //  printf("det=%g\n",temp1->det);  
  //print_matrix(temp1);

  storage[1]=temp1->det;
  

}



void det_remove_anti_fast (MTYPE *storage, double *ts, double *te, int nt, int remid,
		      Matrix * m, Matrix * temp1,Matrix * temp2, Par par){
  storage[0] = m->det;
  //construct space for A^(N+1) and inverse
  int i,j;
  int nsz=nt-1;
  //int sz=ARRAY_SZ;
  int sz=nt;
  
  //  MTYPE * g_tau=par.g_tau;
  //  double beta=par.beta;

  if(nsz>0){  

    //remove column
    //size_t array_sz=sizeof(MTYPE) *sz;
    
    //assert(nsz<ARRAY_SZ);
    MTYPE u[ARRAY_SZ] __attribute__((aligned(64)));
    MTYPE v[ARRAY_SZ] __attribute__((aligned(64)));
    
    //MTYPE *u=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *v=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *u=malloc(sizeof(MTYPE) *sz);
    //MTYPE *v=malloc(sizeof(MTYPE) *sz);
    double tau;
    //MTYPE result;
    
    MTYPE *nm=temp1->mat;
    MTYPE *ni=temp1->inv;
  
    //temp matrix
    MTYPE *am=temp2->mat;
    MTYPE *ai=temp2->inv;

    //original matrix
    MTYPE *mm=m->mat;
    MTYPE *mi=m->inv;

    __assume_aligned(mm,64);
    __assume_aligned(mi,64);
    __assume_aligned(am,64);
    __assume_aligned(ai,64);
    __assume_aligned(nm,64);
    __assume_aligned(ni,64);
    

    int colid = remid;
    int rowid = (remid+1) % sz;

    //remove the id'th column  
#pragma simd
    for (i=0;i<nt;i++){
      u[i]=-mm[i*sz+colid];;
      v[i]=0;
    }
 
    u[rowid]=1.0-mm[rowid*sz+colid];
    v[colid]=1;
    
    det_update_fast(m, temp1, u, v);


    temp1->N=sz;
    temp1->max_sz=ARRAY_SZ;
    
    
    //remove the row
#pragma simd
    for (i=0;i<nt;i++){
      v[i]=-nm[(rowid)*sz+i];
      u[i]=0;
    }

    u[rowid]=1;    
    v[colid]=0;

 
    det_update_fast(temp1, temp2, u, v);
 
    temp1->det=temp2->det;

    temp1->N=nsz;
  
    for(i=0;i<rowid;i++){
#pragma simd
      for(j=0;j<colid;j++){
	nm[i*nsz+j]=am[i*sz+j];
	ni[j*nsz+i]=ai[j*sz+i];
      }
#pragma simd
      for(j=colid;j<nsz;j++){
	nm[i*nsz+j]=am[i*sz+j+1];
	ni[j*nsz+i]=ai[(j+1)*sz+i];
      }
    }
    
    for(i=rowid;i<nsz;i++){
#pragma simd
      for(j=0;j<colid;j++){
	nm[i*nsz+j]=am[(i+1)*sz+j];
	ni[j*nsz+i]=ai[(j)*sz+(i+1)];
      }
#pragma simd
      for(j=colid;j<nsz;j++){
	nm[i*nsz+j]=am[(i+1)*sz+j+1];
	ni[j*nsz+i]=ai[(j+1)*sz+i+1];
      }
    }

    if(rowid==0){
      //move the first column to last
      for(i=0;i<nsz;i++)
	u[i]=nm[i*nsz];//store the first column using u;
      for(i=0;i<nsz;i++)
	for(j=0;j<nsz-1;j++)
	  nm[i*nsz+j]=nm[i*nsz+j+1];
      for(i=0;i<nsz;i++)
	nm[i*nsz+nsz-1]=u[i];

      for(i=0;i<nsz;i++)
	u[i]=ni[i];//store the first row using u;
      for(i=0;i<nsz-1;i++)
	for(j=0;j<nsz;j++)
	  ni[i*nsz+j]=ni[(i+1)*nsz+j];
      for(i=0;i<nsz;i++)
	ni[(nsz-1)*nsz+i]=u[i];


    }

    storage[1]=temp1->det;
  }    
  
  else{
    //new matrix is empty

    temp1->max_sz=ARRAY_SZ;
    temp1->N=0;
    temp1->det = 1.0;
    storage[1]=1.0;
  }

}





void back_det_remove_anti_fast (MTYPE *storage, double *ts, double *te, int nt, int remid,
		      Matrix * m, Matrix * temp1,Matrix * temp2, Par par){
  storage[0] = m->det;
  //construct space for A^(N+1) and inverse
  int i,j;
  int nsz=nt-1;

  if (nsz == 0 ){
    temp1->max_sz=ARRAY_SZ;
    temp1->N=0;
    temp1->det = 1.0;
    storage[1]=1.0;
    return;
  }

  //int sz=ARRAY_SZ;
  int sz=nt;
  double beta=par.beta;
  
  //  MTYPE * g_tau=par.g_tau;
  //  double beta=par.beta;

    //remove column
    //size_t array_sz=sizeof(MTYPE) *sz;
    
    //assert(nsz<ARRAY_SZ);
  MTYPE u[ARRAY_SZ] __attribute__((aligned(64)));
  MTYPE v[ARRAY_SZ] __attribute__((aligned(64)));
    
    //MTYPE *u=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *v=(MTYPE *)_mm_malloc(array_sz,64 );
    //MTYPE *u=malloc(sizeof(MTYPE) *sz);
    //MTYPE *v=malloc(sizeof(MTYPE) *sz);
  double tau;
    //MTYPE result;
    
  MTYPE *nm=temp1->mat;
  MTYPE *ni=temp1->inv;
  
    //temp matrix
  MTYPE *a=temp2->mat;
  MTYPE *ai=temp2->inv;

    //original matrix
  MTYPE *mm=m->mat;
  MTYPE *mi=m->inv;

  __assume_aligned(mm,64);
  __assume_aligned(mi,64);
  __assume_aligned(a,64);
  __assume_aligned(ai,64);
  __assume_aligned(nm,64);
  __assume_aligned(ni,64);
    
  double t_end=te[remid];

  
#pragma simd
  for (i=0;i<nt;i++){
      //tau=te[0][remid]-ts[0][i];
      //get_g_tau (g_tau, tau, beta, &result);
      //result=tmp1->mat[i*sz+remid];
    u[i]=-mm[i*sz+remid];;
    v[i]=0;
  }
 
    //tau=te[0][remid]-ts[0][remid];
    //get_g_tau (g_tau, tau, beta, &result);
    //result=tmp1->mat[remid*sz+remid];
  u[remid]=1-mm[remid*sz+remid];
  v[remid]=1;
    
  det_update_fast(m, temp2, u, v);


  temp2->N=sz;
  temp2->max_sz=ARRAY_SZ;
    
    
    //remove the row
#pragma simd
  for (i=0;i<nt;i++){
      //tau=te[0][i]-ts[0][remid];
      //get_g_tau (g_tau, tau, beta, &result);
      //result=tmp1->mat[remid*sz+i];
    v[i]=-a[remid*sz+i];
    u[i]=0;
  }
    
  v[remid]=0;
  u[remid]=1;
 
  det_update_fast(temp2, temp1, u, v);
 
  temp2->det=temp1->det;
    
  temp2->N=nsz;
  
  for(i=0;i<remid;i++){
#pragma simd
    for(j=0;j<remid;j++){
      a[i*nsz+j]=nm[i*sz+j];
      ai[i*nsz+j]=ni[i*sz+j];
    }
#pragma simd
    for(j=remid;j<nsz;j++){
      a[i*nsz+j]=nm[i*sz+j+1];
      ai[i*nsz+j]=ni[i*sz+j+1];
    }
  }
    
  for(i=remid;i<nsz;i++){
#pragma simd
    for(j=0;j<remid;j++){
      a[i*nsz+j]=nm[(i+1)*sz+j];
      ai[i*nsz+j]=ni[(i+1)*sz+j];
    }
#pragma simd
    for(j=remid;j<nsz;j++){
      a[i*nsz+j]=nm[(i+1)*sz+j+1];
      ai[i*nsz+j]=ni[(i+1)*sz+j+1];
    }
  }



  int modid=(remid-1)%(nsz);
#pragma simd
  for (i=0;i<remid;i++){
    double tau=t_end-ts[i];
    MTYPE result;
    get_g_tau (par.g_tau, tau, beta, &result);
    u[i]=result-a[i*nsz+modid];
    v[i]=0;
    }

  for (i=remid;i<nsz;i++){
    double tau=t_end-ts[i+1];
    MTYPE result;
    get_g_tau (par.g_tau, tau, beta, &result);
    u[i]=result-a[i*nsz+modid];
    v[i]=0;
    }

 
  v[modid]=1;
  
  det_update_fast(temp2, temp1, u, v);

  storage[1]=temp1->det;

}

void construct_matrix(double *ts, double *te, int nt, Matrix *temp, Par par){
  int i,j;
  MTYPE result;
  double beta,tau;
  beta = par.beta;
  for (i=0;i<nt;i++){
    for(j=0;j<nt;j++){
      tau=te[j]-ts[i];
      get_g_tau (par.g_tau, tau, beta, &result);
      temp->mat[i*nt+j]=result;
    }
  }
  temp->N=nt;
}
