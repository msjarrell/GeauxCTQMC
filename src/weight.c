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


//extern unsigned int seed[TASK]; 


//implementations

void weight_shift(Mempool *mempool, int i_o, int modid, double news, double newe, Par par)
{

  double time_start = 0.0;
  double time_end = par.beta;
  double ef = par.ef;
  double U = par.u;
  double beta=par.beta;
  
  // variable calculations
  double *ts,*tsb;
  double *te,*teb;
  int nt;
  int ntb;
  Matrix *m;

  int i_o_s;
  if(i_o % 2)
    i_o_s= i_o-1;
  else 
    i_o_s= i_o+1;
  
  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  tsb= mempool->is[i_o_s];
  teb= mempool->ie[i_o_s];
  nt = mempool->n[i_o];
  ntb= mempool->n[i_o_s];
  m=mempool->m[i_o];

  MTYPE storage[2];

  det_shift_fast (storage, ts,te, nt, modid, news, newe, m, mempool->m_tmp1,mempool->m_tmp2,par);

  MTYPE detold = storage[0];
  MTYPE detnew = storage[1];

  double lnew = newe - news;
  if (lnew < 0)
    lnew += beta;

  double lold = te[modid] - ts[modid];
  if (lold < 0)
    lold += beta;


  double wacc = 1.0;
  wacc *= cabs(detnew / detold); 
  double Q=1.0;


  
  if (ntb>0){

    double deltalov = 0;
    double lov1 =   overlap_total ( tsb, teb, ntb, news, newe, beta); 
    double lov2 =   overlap_total ( tsb, teb, ntb, ts[modid], te[modid], beta);  

    deltalov = lov1 - lov2;

    //printf("deltaov=%g\n",deltalov);
    Q = exp (-1.0 * ef * (lnew - lold)) * exp (-1.0 * U * deltalov);

#ifdef DEBUG_VERBOSE
  printf(" Q=%lf\n",Q);
#endif

  }

  else{
    Q = exp ( -1.0 * ef * (lnew - lold)) * (1.0 + exp( -1.0 * beta * ef) * exp( -1.0 * lnew * U)) / (1.0 + exp( -1.0 * beta * ef) * exp( -1.0 * lold * U));
#ifdef DEBUG_VERBOSE
  printf("!WARNING Q=%lf\n",Q);
#endif

  }
  //printf("Q=%g\n",Q);
  
  wacc *= Q;
  //wacc *= SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);

  //printf("Modify ACC weight =  %f\n",wacc);
  //printf ("modify with rate:%lf\n",wacc);
  if ( myrand(&(mempool->seed)) < wacc) {
    //printf("accept modify\n");
    accept_shift (ts, te, modid, news, newe);
    mempool->shift_accpt++;
    Matrix *temp;

    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;
    
  }

}




void weight_add(Mempool *mempool, int i_o, int newid, double news, double newe, Par par)
{
  //  printf("weight_add\n");

  double time_start = 0.0;
  double time_end = par.beta;
  double ef = par.ef;
  double U = par.u;
  double beta=par.beta;

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
  
  //printf("\nio=%d\tios=%d\n",i_o,i_o_s);
  
  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  nt = mempool->n+i_o;

  tsb= mempool->is[i_o_s];
  teb= mempool->ie[i_o_s];
  ntb= mempool->n[i_o_s];

  m=mempool->m[i_o];

  double lmax = 0;

  MTYPE storage[2];

  det_add_fast (storage, ts, te, *nt, newid, news, newe, m, mempool->m_tmp1, mempool->m_tmp2 , par);

  //print_matrix(m);
  //print_matrix(mempool->m_tmp1);

  int i;

  //for(i=0;i<*nt;i++)
  //    printf("%g\t%g\n",ts[i],te[i]);

  if (*nt > 0) {
    int nextid = (newid ) % (*nt);
    
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
  double lnew = newe - news;
  if (lnew < 0)
    lnew += beta;

  double deltalov =  overlap_total( tsb, teb, ntb, news, newe, beta); 
  /*
  double deltalov2 = 0;
  overlap_single( &deltalov2, tsb, teb, ntb, news, newe, beta); 
  if(abs(deltalov-deltalov2)> 1E-6){
    printf("news=%g\tnewe=%g\n",news,newe);
    printf("ov1=%g\tov2=%g\n",deltalov,deltalov2);
    int i;
    for(i=0;i<ntb;i++)
      printf("s=%g\te=%g\n",tsb[i],teb[i]);
  }
  */
  //printf("deltalov=%g",deltalov);
  /*    
  printf("propose add on spin channel %d:\n",i_o);
  printf("news=%lf,newe=%lf\n",news,newe);
  printf("beta=%lf,lmax=%lf,deltalov=%lf\n",beta,lmax,deltalov);
  printf("detold=%g,detnew=%g\n",creal(detold),creal(detnew));
  printf("k=%d\n",k);
  */   

  double wacc = 1.0;		  
  wacc *= ((beta* lmax) / (k + 1.0));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);
  //  printf("detnew=%g\tdetold=%g\n",creal(detnew),creal(detold));
  //printf("wacc=%g\n",wacc);
  wacc *= cabs(detnew / detold);
  //printf("wacc=%g\n",wacc);

  //   wacc*=v2;

  double Q=1.0;
  if(*nt>0){
    if(ntb==0){
      //nt>0,nd=0
      double lold=length_total ( ts, te, *(nt), beta);

      //printf("10\n");
      //printf("lold=%g\nlnew=%g\n",lold,lnew);
      Q = exp(- 1.0* lnew * ef) *(1.0+ exp(-1.0* beta * ef) * exp(-1.0* (lnew + lold)* U))/ (1.0+exp(-1.0* beta * ef)* exp(-1.0* lold * U));
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

      mempool->add_10 ++;
      //Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
    }
    else
      {
	//printf("11\n");
	//printf("lnew=%g\nlov=%g\n",lnew,deltalov);
	//printf("%g\t%g\t%g\t%g\n",ef,lnew,U,deltalov);
	Q = exp (-1.0 * ef * (lnew))  * exp (-1.0 * U * deltalov);
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

	mempool->add_11 ++;
	//Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
      }
  }
  else{
    if(ntb==0){
      //nt=0,nd=0
      //printf("00\n");
      Q = exp(-1.0* lnew * ef) * ( 1.0 + exp(-1.0 *beta * ef) * exp(-1.0* lnew * U)) / (1.0+ 2.0 * exp(-1.0 * beta * ef) + exp(-1.0* beta * (2.0 * ef + U)));
 
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

      mempool->add_00 ++;
    //Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
    }
    else      {
      //nt=0,nd!=0
      double lold=length_total ( tsb, teb, ntb, beta);

      //printf("01\n");
      //printf("overlap=%g\n",deltalov);
      //printf("lold=%g\nlnew=%g\n",lold,lnew);
      Q = exp(- 1.0* lnew * ef) *exp( - 1.0*deltalov * U) / (1.0+ exp(-1.0 *beta * ef) * exp (-1.0* lold *U)); //l_{\bar{sigma}}
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

      mempool->add_01 ++;
      //Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
    }
  }

  //  printf("Q=%g\n",Q);

  wacc *=Q;
  //printf("wacc=%g\n\n",wacc);

  //printf("test %f %g %g \n",wacc,detnew,detold);
  //  if (wacc >= 1.0)
  //  wacc = 1.0;
  //printf("wacc=%lf\n",wacc);

  if ( myrand(&(mempool->seed)) < wacc) {
    assert((*nt+1)<ARRAY_SZ);
    accept_add (ts,te, nt, newid, news, newe);
    mempool->add_accpt++;

    Matrix *temp;
    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;
    
  }
  //  printf("weight_add done\n");
}


/*
void weight_add_anti(Mempool *mempool, int i_o, int newid, double news, double newe, Par par)
{
  //  printf("weight_add_anti\n");


  double time_start = 0.0;
  double time_end = par.beta;
  double ef = par.ef;
  double U = par.u;
  double beta=par.beta;

  double *ts,*tsb;
  double *te,*teb;
  int *nt;
  int ntb;
  Matrix *m;

  
  if (up) {
    ts = mempool->ius;
    te = mempool->iue;
    tsb= mempool -> ids;
    teb= mempool -> ide;
    nt = &(mempool->nu);
    ntb= mempool->nd;
    m=mempool->m_up;
  }
  else {
    ts = mempool->ids;
    te = mempool->ide;
    tsb= mempool -> ius;
    teb= mempool -> iue;
    nt = &(mempool->nd);
    ntb= mempool->nu;
    m=mempool->m_down;
  }


  double lmax = 0;

  MTYPE storage[2];

  det_add_anti_fast (storage, ts, te, *nt, newid, news, newe, m, mempool->m_tmp1, mempool->m_tmp2 , par);

  //print_matrix(m);
  //print_matrix(mempool->m_tmp1);


  lmax=te[newid]-news;
  if(lmax<0) lmax+=beta;
  
  int k = *nt;
  
  MTYPE detnew = storage[1];
  MTYPE detold = storage[0];
  double lnew = newe - news;
  if (lnew < 0)
    lnew += beta;


  double deltalov = 0;

  overlap_single (&deltalov, tsb, teb, ntb, news, newe, beta);  

  
  //printf("propose add on spin channel %d:\n",up);
  //printf("news=%lf,newe=%lf\n",news,newe);
  //printf("beta=%lf,lmax=%lf,deltalov=%lf\n",beta,lmax,deltalov);
  //printf("detold=%lf,detnew=%lf\n",detold,detnew);
  //printf("k=%d\n",k);
  

  double wacc = 1.0;		  
  wacc *= ((beta * lmax) / (k + 1.0));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);
  //  printf("detnew=%g\tdetold=%g\n",creal(detnew),creal(detold));
  //printf("wacc=%g\n",wacc);
  wacc *= cabs(detnew / detold);
  //printf("wacc=%g\n",wacc);

  //   wacc*=v2;


  double Q=1.0;
  if(*nt>0){
    if(ntb==0){
      //nt>0,nd=0
      double lold=0;
      overlap_single (&lold, ts, te, *(nt), 0, beta, beta);
      //printf("lold=%g\nlnew=%g\n",lold,lnew);
      Q = exp(- 1.0* lnew * ef) *(1.0+ exp(-1.0* beta * ef) * exp(-1.0* (lnew + lold)* U))/ (1.0+exp(-1.0* beta * ef)* exp(-1.0* lold * U));
    }
    else
      {
	//printf("lnew=%g\nlov=%g\n",lnew,deltalov);
	Q = exp (-1.0 * ef * (lnew))  * exp (-1.0 * U * deltalov);      
      }
  }
  else{
    if(ntb==0)
      //nt=0,nd=0
      Q = exp(-1.0* lnew * ef) * ( 1.0 + exp(-1.0 *beta * ef) * exp(-1.0* lnew * U)) / (1.0+ 2.0 * exp(-1.0 * beta * ef) + exp(-1.0* beta * (2.0 * ef + U)));
    else      {
      //nt=0,nd!=0
      double lold=0;
      overlap_single (&lold, tsb, teb, ntb, 0, beta, beta);
      //      printf("lold=%g\nlnew=%g\n",lold,lnew);
      Q = exp(- 1.0* lnew * ef) *exp( - 1.0*deltalov * U) / (1.0+ exp(-1.0 *beta * ef) * exp (-1.0* lold *U)); //l_{\bar{sigma}}
    }
  }

  //printf("Q=%g\n",Q);

  wacc *=Q;
  //printf("wacc=%g\n\n",wacc);

    //printf("test %f %g %g \n",wacc,detnew,detold);
  if (wacc >= 1.0)
    wacc = 1.0;
  //printf("wacc=%lf\n",wacc);

  if ( myrand(&(mempool->seed)) < wacc) {
    assert((*nt+1)<ARRAY_SZ);
    accept_add_anti (ts,te, nt, newid, news, newe);
    mempool->add_anti_accpt++;
    Matrix *temp;
    if (up) {
      //printf ("accept add on up channel\n");
      temp=mempool->m_up;
      mempool->m_up=mempool->m_tmp1;
      mempool->m_tmp1=temp;
    }
    else {
      //printf ("accept add on down channel\n");
      temp=mempool->m_down;
      mempool->m_down=mempool->m_tmp1;
      mempool->m_tmp1=temp;
    }
  }
  //  printf("weight_add done\n");
}
*/


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

  double l = te[remid] - ts[remid];
  if (l < 0)
    l += beta;
    //double lold = te[remid] - ts[remid];
  
  double deltalov = overlap_total ( tsb, teb, ntb, ts[remid], te[remid],beta);

  double Q=1.0;

  if(*nt>1){
    int nextid = (remid + 1) % (*nt);
    lmax = ts[nextid] - ts[remid];
    if (lmax < 0)
      lmax += beta;
  }
  else 
    lmax = beta;

  //printf("deltalov=%f\n",deltalov);
  if((*nt) >1){

    if(ntb==0){
      //nt>0,nd=0
      //      double lold= overlap_total ( ts, te, (*nt), 0, beta,beta);      
      double lold= length_total ( ts, te, (*nt),beta);      

      Q=exp(- l * ef) *(1.0+ exp(-beta * ef) * exp(-(lold)* U))/ (1.0+exp(-beta * ef)* exp(-(lold-l) * U));
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

      mempool->rem_10 ++;
      //Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
    }
    else{
      Q= exp (-1.0 * ef * (l))  * exp ( -1.0* U * deltalov);
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

      mempool->rem_11 ++;
    //Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
    }
  }
  else{
    if(ntb==0){
      //nt=1,nd=0
      Q = exp(-1.0* l * ef) * ( 1.0 + exp(-1.0 *beta * ef) * exp(-1.0* l * U)) / (1.0+ 2.0 * exp(-1.0 * beta * ef) + exp(-1.0* beta * (2.0 * ef + U)));
      //Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

      mempool->rem_00 ++;
    }
    else{
      //nt=1,nd>0
      double lold= length_total ( tsb, teb, ntb, beta);
      //printf("overlap=%g\n",deltalov);
      //printf("lold=%g\nlnew=%g\n",lold,lnew);
      Q = exp(- 1.0* l * ef) *exp( - 1.0 * deltalov * U) / (1.0+ exp(-1.0 *beta * ef) * exp (-1.0* lold *U)); //l_{\bar{sigma}}
#ifdef DEBUG_VERBOSE
      printf("Q=%lf\n",Q);
#endif

      mempool->rem_01 ++;
      //Q*=SDL_WEIGHT(mempool->sdl_m,mempool->is,mempool->ie,mempool->n);
    }

  }
  /*          
       printf("propose remove on spin channel %d:\n",i_o);
       printf("rems=%lf,reme=%lf\n",ts[remid],te[remid]);
       printf("beta=%lf,lmax=%lf,deltalov=%lf\n",beta,lmax,deltalov);
       printf("detold=%g,detnew=%g\n",creal(detold),creal(detnew));
       printf("k=%d\n",k);
       printf("Q=%g\n",Q);
  */
    wacc = 1.0;
    
    wacc *= ((k) / (beta * lmax));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);

    wacc *= cabs(detnew / detold);
    
    wacc *= 1.0/Q;
    
    //printf("wacc=%g\n",wacc);
   
  if (myrand(&(mempool->seed)) < wacc) {
    accept_remove (ts, te, nt, remid);
    mempool->remove_accpt++;
    Matrix * temp;
    temp=mempool->m[i_o];
    mempool->m[i_o]=mempool->m_tmp1;
    mempool->m_tmp1=temp;

  }
  //  printf("weight_remove done\n");
}

/*

void
weight_remove_anti(Mempool *mempool, int up, int remid, Par par)
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
  if (up) {
    ts = mempool->ius;
    te = mempool->iue;
    tsb= mempool -> ids;
    teb= mempool -> ide;
    nt = &(mempool->nu);
    ntb = mempool-> nd;
    m = mempool-> m_up;
  }

  else {
    ts = mempool->ids;
    te = mempool->ide;
    tsb= mempool -> ius;
    teb= mempool -> iue;
    nt = &(mempool->nd);
    ntb = mempool-> nu;
    m = mempool-> m_down;
  }

  
  double lmax = 0;
  
  MTYPE storage[2];
  
  det_remove_anti_fast (storage, ts, te, *nt, remid, m , mempool->m_tmp1,mempool->m_tmp2 , par);

  int k = *nt;
 
  MTYPE detnew = storage[1];
  MTYPE detold = storage[0];

  double l = te[remid] - ts[remid];
  if (l < 0)
    l += beta;
    //double lold = te[remid] - ts[remid];
  
  double deltalov = 0;

  overlap_single (&deltalov, tsb, teb, ntb, ts[remid], te[remid],beta);

  deltalov = deltalov;
  
  
  double Q=1.0;


  if(*nt>1){
    int pvs_id = (remid - 1) % (*nt);
    lmax = ts[remid] - te[pvs_id];
    if (lmax < 0)
      lmax += beta;
  }
  else 
    lmax = beta;


  //printf("deltalov=%f\n",deltalov);
  if(ntb==0){
    //nt>0,nd=0
    double lold=0;
    overlap_single (&lold, ts, te, (*nt), 0, beta,beta);      
    Q=exp(- l * ef) *(1.0+ exp(-beta * ef) * exp(-(lold)* U))/ (1.0+exp(-beta * ef)* exp(-(lold-l) * U));
  }
  else
    Q= exp (-1.0 * ef * (l))  * exp ( -1.0* U * deltalov);      
  

    wacc= 1.0;
    
    //        wacc/=v2;

    wacc *= ((k) / (beta * lmax));	//printf("test %f %i %f %i \n",wacc,k,lmax,nt);
    wacc *= cabs(detnew / detold);
    //printf("test %f %g %g \n",wacc,detnew,detold);
    wacc /= Q;

    //    wacc /=v2;
    //printf("test %f %g %g \n",wacc,detnew,detold);
    if (wacc >= 1)
      wacc = 1.0;
    //printf("wacc=%lf\n",wacc);    

    
    

  if (myrand(&(mempool->seed)) < wacc) {
    accept_remove_anti (ts, te, nt, remid);
    mempool->rem_anti_accpt++;
    Matrix * temp;
    if (up) {
      temp=mempool->m_up;
      mempool->m_up=mempool->m_tmp1;
      mempool->m_tmp1=temp;
    }
    else {
      temp=mempool->m_down;
      mempool->m_down=mempool->m_tmp1;
      mempool->m_tmp1=temp;

    }
  }
  //  printf("weight_remove done\n");
}
*/
