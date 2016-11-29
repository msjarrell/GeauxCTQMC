#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <limits.h>
#include <assert.h>

#include "matrix.h"
#include "utility.h"
#include "accept.h"
#include "main.h"
#include "det.h"
#include "parameters.h"

#ifdef SDF
#include <SDF.h>
#endif

//extern unsigned int seed[TASK]; 




void init_hamiltonian( Par *par){
  //  printf("Initializing the hamiltonian\n");

  FILE *fp=fopen("hamiltonian_aim.txt","r");
  //FILE *fp=fopen("hamiltonian_dhm.txt","r");
  int i,j;
  int n     = HM_SZ;
  double hm[HM_SZ*HM_SZ];
  for(i=0;i<HM_SZ;i++){
    for(j=0;j<HM_SZ;j++){
      fscanf(fp,"%lf\t", hm+i*HM_SZ+j);
    }
  }
  fclose(fp);


#ifdef DEBUG_VERBOSE
  printf("Hamiltonian:\n");
  for(i=0;i<HM_SZ;i++){
    for(j=0;j<HM_SZ;j++)
      printf("%lf\t",hm[i*n+j]);
    printf("\n");
  }
#endif



  int lda   = HM_SZ;
  char jobz  =  'V';
  char uplo  =  'U';
  double *w     = (double *)malloc(sizeof(double)*n);
  LAPACKE_dsyev(LAPACK_ROW_MAJOR,jobz, uplo, n, hm, lda, w);


  for(i=0;i<n;i++){
    par->V[i]=w[i];
  }

#ifdef DEBUG_VERBOSE
  printf("Eigen vector:\n");
  for(i=0;i<n;i++)
      printf("%lf\n",w[i]);
#endif


  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      par->U[i*n+j]=hm[i*n+j];
    }

#ifdef DEBUG_VERBOSE
  printf("Eigen matrix:\n");
  for(i=0;i<HM_SZ;i++){
    for(j=0;j<HM_SZ;j++)
      printf("%lf\t",hm[i*n+j]);
    printf("\n");
  }
#endif


  int cd_up_matrix[HM_SZ][HM_SZ];
  int cd_dn_matrix[HM_SZ][HM_SZ];

  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      cd_up_matrix[i][j]=0;
      cd_dn_matrix[i][j]=0;
    }

  /* operators for Anderson Model  */  
  cd_up_matrix[1][0] = 1;
  cd_up_matrix[3][2] = 1;
    
  cd_dn_matrix[2][0] = 1;
  cd_dn_matrix[3][1] = 1;

  /* operators for DHM*/
  /*
  cd_up_matrix[5][4] = 1;
  cd_up_matrix[7][6] = 1;
    
  cd_dn_matrix[6][4] = 1;
  cd_dn_matrix[7][5] = 1;
  */
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      par->OP[0][i*n+j]=cd_up_matrix[i][j];

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      par->OP[1][i*n+j]=cd_up_matrix[j][i];

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      par->OP[2][i*n+j]=cd_dn_matrix[i][j];

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      par->OP[3][i*n+j]=cd_dn_matrix[j][i];

  
  double temp[HM_SZ*HM_SZ];
  
  int m=HM_SZ;

  for(i=0;i<4;i++){
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,  m, m, m, 1.0, par->OP[i], m, hm, m, 0.0, temp, m);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,  m, m, m, 1.0, hm, m, temp, m, 0.0, par->OP[i], m);
  }
  
  free(w);  


#ifdef SDF

  char* pEnd;
  double* vecd=new double[HM_SZ]();
  char mystring[100];
  const char* const file_name = "./SDF_exp_API_native/vecd.txt";

  FILE *f = fopen(file_name, "r");
  int k=0;
  while( fgets(mystring, sizeof(mystring), f) ){
    par->vecd[k] = strtod(mystring, &pEnd);
    k++;
  }
  fclose(f);

#endif


}


void exponential_matrix( double *v, double t,double *expvt){
  // expvt = exp(v * t)
  int i;
  int m=HM_SZ;
#pragma simd
  for(i=0;i<m;i++){
    expvt[i*m+i]=exp(-v[i]*t);
  }
}


void copy_matrix( double *m1, double *m2){
  // m2=m1
  int i;
#pragma simd
  for(i=0;i< (HM_SZ * HM_SZ );i++)
    m2[i]=m1[i];
}

void mmmul( double *m1, double *m2){
  // m2 = m1 * m2
  double temp[HM_SZ*HM_SZ];
  int m=HM_SZ;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,  m, m, m, 1.0, m1, m, m2, m, 0.0, temp, m);
  copy_matrix(temp,m2);
}


void mvmul( double *m1, double *v, double *v_out){
  //v_out = m1 * v 
  cblas_dgemv(CblasRowMajor, CblasNoTrans, HM_SZ, HM_SZ, 1.0, m1, HM_SZ, v, 1, 0.0, v_out, 1);  
}

#ifdef SDF
double weight_loc_SDF(int n, int *op_list,double *time_list, Par par){
  //using krylov method to calculate the trace
  double v2[HM_SZ];
  double inVec[HM_SZ];
  double vecd[HM_SZ];
  double trace = 0;
  int i;
  int i_eig,i_op;
  int KRYLOV_SZ=3;
  double beta = par.beta;

  for( i = 0; i < HM_SZ; i++)
    vecd[i] = par.vecd[i];

  for(i_eig = 0; i_eig < HM_SZ; i_eig ++){
    for( i = 0; i < HM_SZ; i++)
      inVec[i] = par.U[i_eig * HM_SZ + i];
    
    double t1=0;
    double t2;

    for(i_op=0;i_op<n;i_op++){
      t2=time_list[i_op];
      mat_exp_krylov(HM_SZ, KRYLOV_SZ, (t2-t1), inVec, vecd, v2);
      mvmul(par.OP[op_list[i_op]], v2, inVec);
      t1 = t2;
    }
  
    mat_exp_krylov(HM_SZ, KRYLOV_SZ, (beta-t1), inVec, vecd, v2);
    
    for( i = 0; i < HM_SZ; i++)
      trace += par.U[i_eig * HM_SZ + i] * v2[i];
    
  }

  return trace;
}
#endif


double weight_loc(int n, int *op_list,double *time_list, Par par){
  //using matrix-matrix multiplication to calculate the trace
  double *v = par.V;
  double *U = par.U;
  double mloc[HM_SZ*HM_SZ];
  double evt[HM_SZ*HM_SZ] = {0};
  double trace=0;
  int i,j;
  int i_op;


  int m = HM_SZ;

  for(i=0;i<m;i++)
#pragma simd
    for(j=0;j<m;j++)
      mloc[i*m+j]=U[j*m+i];
  
#ifdef DEBUG_VERBOSE_TRACE
  printf("---------------\n");
  for(i=0;i<HM_SZ;i++){
    for(j=0;j<HM_SZ;j++)
      printf("%lf\t",mloc[i*m+j]);
    printf("\n");
  }
#endif

  double t1=0;
  double t2;


  for(i_op=0;i_op<n;i_op++){
    t2=time_list[i_op];
  
    exponential_matrix(v,(t2-t1),evt);
    mmmul( evt, mloc);


#ifdef DEBUG_VERBOSE_TRACE
    printf("t2=%lf\n",t2);
    printf("---------------\n");
    for(i=0;i<HM_SZ;i++){
      for(j=0;j<HM_SZ;j++)
	printf("%lf\t",mloc[i*m+j]);
      printf("\n");
    }
#endif

    mmmul( par.OP[op_list[i_op]],mloc);

    t1=t2;
    
#ifdef DEBUG_VERBOSE_TRACE
    printf("---------------\n");
    for(i=0;i<HM_SZ;i++){
      for(j=0;j<HM_SZ;j++)
	printf("%lf\t",mloc[i*m+j]);
      printf("\n");
    }
#endif
    
  }

  t2=par.beta;
  exponential_matrix(v,(t2-t1),evt);
  mmmul( evt, mloc);
  mmmul( U, mloc);


#ifdef DEBUG_VERBOSE_TRACE
  printf("---------------\n");
  for(i=0;i<HM_SZ;i++){
    for(j=0;j<HM_SZ;j++)
      printf("%lf\t",mloc[i*m+j]);
    printf("\n");
  }
#endif

  
  for(i=0;i<m;i++)
    trace += mloc[i*m+i];


#ifdef DEBUG_VERBOSE_TRACE
  printf("trace=%lf\n",trace);
#endif
  return trace;
}

void generate_rem_list(int n,int *op_list, int *tmp_op_list,double *time_list, double * tmp_time_list, int i_o, double olds, double olde){

#ifdef DEBUG_VERBOSE
    printf("olds=%lf,olde=%lf\n",olds,olde);
#endif


  if(n==2)
    return;
  
  int i;
  int op1,op2;
  double time1,time2;
    
  if(olds<olde){
    op1 = 2* i_o+1;
    op2 = 2*i_o;
    time1 = olde;
    time2 = olds;
  }
  else{
    op1= 2*i_o;
    op2= 2*i_o+1;
    time1 = olds;
    time2 = olde;
  }
  

  i = n-1;
  while((time_list[i] != time1) || (op_list[i] != op1)){
    tmp_time_list[i-2] = time_list[i];
    tmp_op_list[i-2] = op_list[i];
    i--;
#ifdef DEBUG_VERBOSE
    if(i<0){
      printf("Failure!\n");
      printf("To remove: %lf\t%lf\nList:\n",olds,olde);
      for(i=0;i<n;i++)
	printf("%d\t%lf\n",op_list[i],time_list[i]);
      exit(0);
    }
#endif
  }

  i--;
  if(i<0)
    return;

  while((time_list[i] !=time2) || (op_list[i] != op2)){
    tmp_time_list[i-1] = time_list[i];
    tmp_op_list[i-1] = op_list[i];
    i--;

#ifdef DEBUG_VERBOSE
    if(i<0){
      printf("Failure!\n");
      printf("To remove: %lf\t%lf\nList:\n",olds,olde);
      for(i=0;i<n;i++)
	printf("%d\t%lf\n",op_list[i],time_list[i]);
      exit(0);
    }
#endif

  }

  i--;
  if(i<0)
    return;

  
  while(i != -1){
    tmp_time_list[i] = time_list[i];
    tmp_op_list[i] = op_list[i];
    i--;
  }

}



void generate_shift_list(int n,int *op_list, int *tmp_op_list,double *time_list, double * tmp_time_list, int op, double oldt, double newt){
  int i=0;
  if (newt>oldt){
    while((time_list[i] != oldt) || (op_list[i] != op )){
      tmp_time_list[i]=time_list[i];
      tmp_op_list[i]=op_list[i];
      i++;
    }
    
    while((time_list[i+1] < newt) && (i!=n-1)){
      tmp_time_list[i]=time_list[i+1];
      tmp_op_list[i]=op_list[i+1];
      i++;
    }

    tmp_time_list[i]=newt;
    tmp_op_list[i]= op;
    i++;

    while(i<n){
      tmp_time_list[i]=time_list[i];
      tmp_op_list[i]=op_list[i];
      i++;
    }

  }
  else{
    while((time_list[i] < newt) ){
      tmp_time_list[i]=time_list[i];
      tmp_op_list[i]=op_list[i];
      i++;
    }

    tmp_time_list[i]=newt;
    tmp_op_list[i]= op;    
    i++;

    while((time_list[i-1] !=oldt) || (op_list[i-1] != op )){
      tmp_time_list[i]=time_list[i-1];
      tmp_op_list[i]=op_list[i-1];
      i++;    

    }
    while(i!=n){
      tmp_time_list[i]=time_list[i];
      tmp_op_list[i]=op_list[i];
      i++;
    }

  }
  
}


void generate_add_list(int n,int *op_list,double *time_list, int *tmp_op_list, double * tmp_time_list, int i_o, double news, double newe){

#ifdef DEBUG_VERBOSE
    printf("news=%lf,newe=%lf\n",news,newe);
#endif

  if(n==0){
    //original list is empty
    if(news<newe){
      tmp_op_list[0]= 2*i_o;
      tmp_op_list[1]= 2*i_o+1;
      tmp_time_list[0]=news;
      tmp_time_list[1]=newe;
    }
    else{
      tmp_op_list[0]= 2*i_o+1;
      tmp_op_list[1]= 2*i_o;
      tmp_time_list[0]=newe;
      tmp_time_list[1]=news;
    }

#ifdef DEBUG_VERBOSE
    printf("added operators:\n");
    printf("op#1=%d\t%lf:\n",tmp_op_list[0],tmp_time_list[0]);
    printf("op#2=%d\t%lf:\n",tmp_op_list[1],tmp_time_list[1]);
#endif
  }
  else{
    //original list is not empty, copy and insert
    int op1,op2;
    double time1,time2;
    
    if(news<newe){
      op1 = 2* i_o+1;
      op2 = 2*i_o;
      time1 = newe;
      time2 = news;
    }
    else{
      op1= 2*i_o;
      op2= 2*i_o+1;
      time1 = news;
      time2 = newe;
    }
    
    int i= n-1;
    while(time_list[i]>time1){
      tmp_time_list[i+2]=time_list[i];
      tmp_op_list[i+2]=op_list[i];
      i--;
    }

    tmp_op_list[i+2]= op1;
    tmp_time_list[i+2] = time1;
    
    while(time_list[i]>time2){
      tmp_time_list[i+1]=time_list[i];
      tmp_op_list[i+1]=op_list[i];
      i--;
    }

    tmp_op_list[i+1]= op2;
    tmp_time_list[i+1] = time2;

    while(i != -1){
      tmp_time_list[i]=time_list[i];
      tmp_op_list[i]=op_list[i];
      i--;
    }

  }  

}

void weight_shift(Mempool *mempool, int i_o, int modid, double news, double newe, Par par)
{

#ifdef DEBUG_VERBOSE_SHIFT
  printf("modid=%d,\tnews=%lf\t,newe=%lf\n",modid,news,newe);
#endif
  double time_start = 0.0;
  double time_end = par.beta;
  double ef = par.ef;
  double U = par.u;
  double beta=par.beta;
  
  // variable calculations
  double *ts;
  double *te;
  int nt;
  Matrix *m;

  int *op_list = mempool->op_list;
  double *time_list = mempool->time_list;
  int *tmp_op_list = mempool->tmp_op_list;
  double *tmp_time_list = mempool->tmp_time_list;

  
  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  nt = mempool->n[i_o];
  m=mempool->m[i_o];

#ifdef DEBUG_SHIFT
  double *tsb;
  double *teb;
  int ntb;

  int i_o_s;
  if(i_o % 2)
    i_o_s= i_o-1;
  else 
    i_o_s= i_o+1;

  tsb= mempool->is[i_o_s];
  teb= mempool->ie[i_o_s];
  ntb= mempool->n[i_o_s];



  double lnew = newe - news;
  if (lnew < 0)
    lnew += beta;

  double lold = te[modid] - ts[modid];
  if (lold < 0)
    lold += beta;

#endif


  MTYPE storage[2];

  det_shift_fast (storage, ts,te, nt, modid, news, newe, m, mempool->m_tmp1,mempool->m_tmp2,par);

  MTYPE detold = storage[0];
  MTYPE detnew = storage[1];

  double wacc = 1.0;
  wacc *= cabs(detnew / detold); 

  double old_trace = mempool->trace;
  
  int n=mempool->n_sum;

  double newt,oldt;
  int op;

  if(ts[modid] == news){
    newt = newe;
    oldt = te[modid];
    op = i_o * 2 + 1;
  }
  else{
    newt = news;
    oldt = ts[modid];
    op = i_o*2;
  }

  generate_shift_list(n,op_list,tmp_op_list,time_list,tmp_time_list, op, oldt, newt);


#ifdef DEBUG_VERBOSE_SHIFT
  printf("old list:\n");
  print_list(n,op_list,time_list);
  printf("new list:\n");
  print_list(n,tmp_op_list,tmp_time_list);
#endif

  
  double new_trace = weight_loc(n,tmp_op_list,tmp_time_list,par);

#ifndef DEBUG_SHIFT
  double Q = (new_trace/old_trace);
#else
  double Q1 = (new_trace/old_trace);

  double Q;

  if (ntb>0){

    double deltalov = 0;
    double lov1 =   overlap_total ( tsb, teb, ntb, news, newe, beta); 
    double lov2 =   overlap_total ( tsb, teb, ntb, ts[modid], te[modid], beta);  

    deltalov = lov1 - lov2;

    //printf("deltaov=%g\n",deltalov);
    Q = exp (-1.0 * ef * (lnew - lold)) * exp (-1.0 * U * deltalov);


  }

  else{
    double l_tot = overlap_total(ts,te,nt,0,beta,beta);
    
    Q = exp ( -1.0 * ef * (lnew - lold)) * (1.0 + exp( -1.0 * beta * ef) * exp( -1.0 * (l_tot-lold+lnew) * U)) / (1.0 + exp( -1.0 * beta * ef) * exp( -1.0 * l_tot * U));
    
    //Q = exp ( -1.0 * ef * (lnew - lold)) * (1.0 + exp( -1.0 * beta * ef) * exp( -1.0 * lnew * U)) / (1.0 + exp( -1.0 * beta * ef) * exp( -1.0 * lold * U));

  }

  if(fabs((Q1 - Q) / Q)> 1e-3 ){
    printf("modid=%d,\tnews=%lf\t,newe=%lf\n",modid,news,newe);
    printf("modid=%d,\tolds=%lf\t,olde=%lf\n",modid,ts[modid],te[modid]);
    printf("Q =%lf\n",Q);
    printf("Q1=%lf\n",Q1);
    printf("new_trace=%lf\n",new_trace);
    printf("old_trace=%lf\n",old_trace);
    printf("old list:\n");
    print_list(n,op_list,time_list);
    printf("new list:\n");
    print_list(n,tmp_op_list,tmp_time_list);
    getchar();
}

#endif

#ifdef DEBUG_VERBOSE
  printf("Q=%lf\n",Q);
#endif

  wacc *= Q;

  //wacc *= SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);

  //printf("Modify ACC weight =  %f\n",wacc);
  //printf ("modify with rate:%lf\n",wacc);
  if ( myrand(&(mempool->seed)) < wacc) {
    //printf("accept modify\n");
    accept_shift (ts, te, modid, news, newe);
    mempool->shift_accpt++;
    
    mempool->trace=new_trace;
    
    Matrix *temp;
    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;

    mempool->tmp_op_list = op_list;
    mempool->op_list = tmp_op_list;
    mempool->time_list = tmp_time_list;
    mempool->tmp_time_list = time_list;


#ifdef DEBUG_VERIFY
    verify_matrix(i_o,  mempool,  par);
#endif

    
  }

}




void weight_add(Mempool *mempool, int i_o, int addid, double news, double newe, Par par)
{
  //  printf("weight_add\n");

  double time_start = 0.0;
  double time_end = par.beta;
  double ef = par.ef;
  double U = par.u;
  double beta=par.beta;

  double *ts;
  double *te;
  int *nt;
  Matrix *m;
    
  int *op_list = mempool->op_list;
  double *time_list = mempool->time_list;
  int *tmp_op_list = mempool->tmp_op_list;
  double *tmp_time_list = mempool->tmp_time_list;



  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  nt = mempool->n+i_o;

  m=mempool->m[i_o];

  double lmax = 0;

  MTYPE storage[2];

  det_add_fast (storage, ts, te, *nt, addid, news, newe, m, mempool->m_tmp1, mempool->m_tmp2 , par);

  //print_matrix(m);
  //print_matrix(mempool->m_tmp1);

  int i;

  //for(i=0;i<*nt;i++)
  //    printf("%g\t%g\n",ts[i],te[i]);

  if (*nt > 0) {
    int nextid = (addid ) % (*nt);
    
    //printf("newid=%d\tnextid=%d\n",newid,nextid);
    //printf("te_current=%g\nts_next=%g\n",te[newid],ts[nextid]);
    
    lmax = ts[nextid] - news;
    if (lmax < 0)
      lmax += beta;
  }
  else {
    lmax=beta;
  }

  //  printf("news=%g\tnewe=%g\n",news,newe);
  //printf("lmax=%g\n",lmax);
  
  int k = *nt;
  
  MTYPE detnew = storage[1];
  MTYPE detold = storage[0];

  double wacc = 1.0;		  
  wacc *= ((beta* lmax) / (k + 1.0));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);
  wacc *= cabs(detnew / detold);

  double old_trace = mempool->trace;
    
  int n = mempool->n_sum;

  generate_add_list(n,op_list,time_list, tmp_op_list, tmp_time_list,i_o,  news, newe); 


#ifdef DEBUG_VERBOSE
  printf("old list:\n");
  print_list(n,op_list,time_list);
  printf("new list:\n");
  print_list(n+2,tmp_op_list,tmp_time_list);
#endif

  double new_trace = weight_loc( (n+2),tmp_op_list,tmp_time_list,par);

  double Q = (new_trace/old_trace);
#ifdef DEBUG_VERBOSE
  printf("Q=%lf\n",Q);
#endif
  wacc *= Q;

 

  if ( myrand(&(mempool->seed)) < wacc) {
    assert((*nt+1)<ARRAY_SZ);

    accept_add (ts,te, nt, addid, news, newe);
    mempool->add_accpt++;
    mempool->trace=new_trace;
    mempool->n_sum += 2;
    /*
    int *tmp1=op_list;
    op_list=tmp_op_list;
    tmp_op_list=tmp1;

    double *tmp2=time_list;
    time_list=tmp_time_list;
    tmp_time_list=tmp2;
    */
    mempool->tmp_op_list = op_list;
    mempool->op_list = tmp_op_list;
    mempool->time_list = tmp_time_list;
    mempool->tmp_time_list = time_list;

    Matrix *temp;
    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;


#ifdef DEBUG_VERIFY
    verify_matrix(i_o,  mempool,  par);
#endif



  }
  //  printf("weight_add done\n");
}


void
weight_remove(Mempool *mempool, int i_o, int remid, Par par)
{

  //printf("weight_remove\n");
  double wacc = 1.0;		//fabs(tau - tauprime);     printf("test %f %f %f \n",wacc,tau,tauprime);

  //temp inputs
  double beta = par.beta; 
  double time_start = 0.0;
  double time_end = beta;
  double ef = par.ef;
  double U = par.u;

  // variable calculations
  double *ts,*tsb;
  double *te,*teb;
  int *nt;
  int ntb;
  Matrix *m;

  int i_o_s;
  if(i_o % 2)
    i_o_s= i_o-1;
  else 
    i_o_s= i_o+1;
  //printf("io=%d\tios=%d\n",i_o,i_o_s);

  int *op_list = mempool->op_list;
  double *time_list = mempool->time_list;
  int *tmp_op_list = mempool->tmp_op_list;
  double *tmp_time_list = mempool->tmp_time_list;

  
  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  tsb= mempool -> is[i_o_s];
  teb= mempool -> ie[i_o_s];
  nt = mempool->n+i_o;
  ntb= mempool->n[i_o_s];
  m=mempool->m[i_o];
  
  double lmax = 0;
  
  MTYPE storage[2];
  
  det_remove_fast (storage, ts, te, *nt, remid, m , mempool->m_tmp1,mempool->m_tmp2 , par);

  int k = *nt;
 
  MTYPE detnew = storage[1];
  MTYPE detold = storage[0];

  if(*nt>1){
    int nextid = (remid + 1) % (*nt);
    lmax = ts[nextid] - ts[remid];
    if (lmax < 0)
      lmax += beta;
  }
  else 
    lmax = beta;

    wacc = 1.0;    
    wacc *= ((k) / (beta * lmax));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);
    wacc *= cabs(detnew / detold);
    
    double old_trace = mempool->trace;
    int n=mempool->n_sum;
    double rems=ts[remid];
    double reme=te[remid];
    generate_rem_list(n,op_list,tmp_op_list,time_list,tmp_time_list, i_o, rems, reme); 

#ifdef DEBUG_VERBOSE
  printf("old list:\n");
  print_list(n,op_list,time_list);
  printf("new list:\n");
  print_list(n-2,tmp_op_list,tmp_time_list);
#endif

 
  double new_trace = weight_loc( (n-2),tmp_op_list,tmp_time_list,par);

  double Q = (new_trace/old_trace);
#ifdef DEBUG_VERBOSE
  printf("Q=%lf\n",1.0/Q);
#endif
  wacc *= Q;

    
   
  if (myrand(&(mempool->seed)) < wacc) {
    accept_remove (ts, te, nt, remid);
    mempool->remove_accpt++;
    mempool->trace=new_trace;
    Matrix * temp;
    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;

    mempool->n_sum -= 2;

    mempool->tmp_op_list = op_list;
    mempool->op_list = tmp_op_list;
    mempool->time_list = tmp_time_list;
    mempool->tmp_time_list = time_list;


#ifdef DEBUG_VERIFY
    verify_matrix(i_o,  mempool,  par);
#endif


  }
  //  printf("weight_remove done\n");
}


void weight_add_anti(Mempool *mempool, int i_o, int addid, double news, double newe, Par par)
{
  //  printf("weight_add\n");

  double time_start = 0.0;
  double time_end = par.beta;
  double ef = par.ef;
  double U = par.u;
  double beta=par.beta;

  double *ts;
  double *te;
  int *nt;
  Matrix *m;
    
  int *op_list = mempool->op_list;
  double *time_list = mempool->time_list;
  int *tmp_op_list = mempool->tmp_op_list;
  double *tmp_time_list = mempool->tmp_time_list;


  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  nt = mempool->n+i_o;

  m=mempool->m[i_o];

  double lmax = 0;

  MTYPE storage[2];

  det_add_anti_fast (storage, ts, te, *nt, addid, news, newe, m, mempool->m_tmp1, mempool->m_tmp2 , par);

  //print_matrix(m);
  //print_matrix(mempool->m_tmp1);

  int i;

  if (*nt > 0) {
    int nextid = addid;
#ifdef DEBUG_VERBOSE
    printf("newid=%d\tnextid=%d\n",addid,nextid);
#endif
    //printf("te_current=%g\nts_next=%g\n",te[newid],ts[nextid]);
    
    lmax = te[nextid] - news;
    if (lmax < 0)
      lmax += beta;
  }
  else {
    lmax=beta;
  }

  //  printf("news=%g\tnewe=%g\n",news,newe);
  //printf("lmax=%g\n",lmax);
  
  int k = *nt;
  
  MTYPE detnew = storage[1];
  MTYPE detold = storage[0];

  double wacc = 1.0;		  
  wacc *= ((beta* lmax) / (k + 1.0));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);
  wacc *= cabs(detnew / detold);

  double old_trace = mempool->trace;
    
  int n = mempool->n_sum;

  generate_add_list(n,op_list,time_list, tmp_op_list, tmp_time_list,i_o,  newe, news); 


#ifdef DEBUG_VERBOSE
  printf("old list:\n");
  print_list(n,op_list,time_list);
  printf("new list:\n");
  print_list(n+2,tmp_op_list,tmp_time_list);
#endif

  double new_trace = weight_loc( (n+2),tmp_op_list,tmp_time_list,par);

  double Q = (new_trace/old_trace);
#ifdef DEBUG_VERBOSE
  printf("Q=%lf\n",Q);
#endif
  wacc *= Q;

 
  if ( myrand(&(mempool->seed)) < wacc) {
    assert((*nt+1)<ARRAY_SZ);
    mempool->isFull[i_o] = 0;
    accept_add_anti (ts,te, nt, addid, news, newe);
    mempool->add_anti_accpt++;
    mempool->trace=new_trace;
    mempool->n_sum += 2;
    /*
    int *tmp1=op_list;
    op_list=tmp_op_list;
    tmp_op_list=tmp1;

    double *tmp2=time_list;
    time_list=tmp_time_list;
    tmp_time_list=tmp2;
    */
    mempool->tmp_op_list = op_list;
    mempool->op_list = tmp_op_list;
    mempool->time_list = tmp_time_list;
    mempool->tmp_time_list = time_list;

    Matrix *temp;
    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;


#ifdef DEBUG_VERIFY
    verify_matrix( i_o, mempool,par);
#endif


  }
  //  printf("weight_add done\n");
}



void
weight_remove_anti(Mempool *mempool, int i_o, int remid, Par par)
{

  //printf("weight_remove\n");
  double wacc = 1.0;		//fabs(tau - tauprime);     printf("test %f %f %f \n",wacc,tau,tauprime);

  //temp inputs
  double beta = par.beta; 
  double time_start = 0.0;
  double time_end = beta;
  double ef = par.ef;
  double U = par.u;

  // variable calculations
  double *ts,*tsb;
  double *te,*teb;
  int *nt;
  int ntb;
  Matrix *m;

  int i_o_s;
  if(i_o % 2)
    i_o_s= i_o-1;
  else 
    i_o_s= i_o+1;
  //printf("io=%d\tios=%d\n",i_o,i_o_s);

  int *op_list = mempool->op_list;
  double *time_list = mempool->time_list;
  int *tmp_op_list = mempool->tmp_op_list;
  double *tmp_time_list = mempool->tmp_time_list;

  
  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  tsb= mempool -> is[i_o_s];
  teb= mempool -> ie[i_o_s];
  nt = mempool->n+i_o;
  ntb= mempool->n[i_o_s];
  m=mempool->m[i_o];
  
  double lmax = 0;
  
  MTYPE storage[2];
  
  det_remove_anti_fast (storage, ts, te, *nt, remid, m , mempool->m_tmp1,mempool->m_tmp2 , par);

  int k = *nt;
 
  MTYPE detnew = storage[1];
  MTYPE detold = storage[0];

  if(*nt>1){
    int nextid = (remid + 1) % (*nt);
    lmax = te[nextid] - te[remid];
    if (lmax < 0)
      lmax += beta;
  }
  else 
    lmax = beta;

    wacc = 1.0;    
    wacc *= ((k) / (beta * lmax));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);
    wacc *= cabs(detnew / detold);
    
    double old_trace = mempool->trace;
    int n=mempool->n_sum;
    double rems=ts[(remid+1)%k];
    double reme=te[remid];
    generate_rem_list(n,op_list,tmp_op_list,time_list,tmp_time_list, i_o, rems, reme); 

#ifdef DEBUG_VERBOSE
    printf("old list:\n");
    print_list(n,op_list,time_list);
    printf("new list:\n");
    print_list(n-2,tmp_op_list,tmp_time_list);
#endif

 
  double new_trace = weight_loc( (n-2),tmp_op_list,tmp_time_list,par);

  double Q = (new_trace/old_trace);
#ifdef DEBUG_VERBOSE
  printf("Q=%lf\n",1.0/Q);
#endif
  wacc *= Q;

    
   
  if (myrand(&(mempool->seed)) < wacc) {
    accept_remove_anti (ts, te, nt, remid);
    mempool->rem_anti_accpt++;
    mempool->trace=new_trace;
    Matrix * temp;
    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;

    mempool->n_sum -= 2;

    mempool->tmp_op_list = op_list;
    mempool->op_list = tmp_op_list;
    mempool->time_list = tmp_time_list;
    mempool->tmp_time_list = time_list;

    if (mempool->n[i_o] == 0)
      mempool->isFull[i_o] = 1;

#ifdef DEBUG_VERIFY
    verify_matrix( i_o, mempool, par);
#endif

  }
  //  printf("weight_remove done\n");
}

