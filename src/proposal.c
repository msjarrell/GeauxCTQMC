#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"
#include "matrix.h"
#include "proposal.h"
#include "utility.h"
#include "weight_dhm.h"


void proposal_next(Mempool *mempool, Par par)
{

  double r = myrand(&(mempool->seed));

  int CHOICES=5;
  //2: add/remove only
  //3: modify/add/remove of segments
  //5: modify/add/remove of segments, add/remove of anti-segments
  int choice= floor(CHOICES * r) ;

  mempool->iter++;

  switch(choice){
 
  case 0:
#ifdef DEBUG_VERBOSE
    printf("Adding...\n");
#endif
    mempool->add_prop++;
    proposal_add( mempool, par);
    break;

  case 1:
#ifdef DEBUG_VERBOSE
    printf("Removing...\n");
#endif
    mempool->remove_prop++;
    proposal_remove(mempool, par);
    break;

  case 2:{
#ifdef DEBUG_VERBOSE
    printf("Shifting...\n");
#endif
    //int i;
    //for(i=0;i<10;i++){
    mempool->shift_prop++;
    proposal_shift(mempool, par);
      //}
  }
    break;

  case 3:
#ifdef DEBUG_VERBOSE
    printf("Adding anti...\n");
#endif
    mempool->add_anti_prop++;
    proposal_add_anti(mempool,par);
    break;

  case 4:
#ifdef DEBUG_VERBOSE
    printf("Removing anti...\n");
#endif
    mempool->rem_anti_prop++;
    proposal_rem_anti(mempool,par);

  }

}


void
proposal_shift (Mempool *mempool, Par par)
{

  double time_start = 0.0;
  double time_end = par.beta;
  double beta=par.beta;

  double r = myrand(&(mempool->seed));
  
  int i_o = 0;
  int modid;
  double news,newe;
  double *ts;
  double *te;
  int nt;
  /*
  int n_sum=0;


  for(i_o=0;i_o<N_ORBIT*2;i_o++)
    n_sum +=  mempool->n[i_o];

  if (n_sum == 0){
#ifdef DEBUG_VERBOSE
    printf("s_sum=0,returning...\n");
#endif

    return;
  }

  modid = floor(r * n_sum);

  i_o = 0;
  nt = mempool->n[i_o];
  while(modid >= nt){
    modid -= nt;
    i_o ++;
    nt = mempool->n[i_o];
  }
  */

  r = myrand(&(mempool->seed));  
  i_o = floor(N_CHANNEL * r);

  r = myrand(&(mempool->seed));  
  
  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  nt = mempool->n[i_o];  

  if (nt==0)
    return;

  modid = floor( nt * r);

  
  double left;
  double right;

  int pvsid = modid - 1;
  if (pvsid<0) pvsid+=nt;
  int nextid = (modid + 1)%(nt);

  //  printf("id=%d\tpvs=%d\tnext=%d\n",modid,pvsid,nextid);

  left = te[pvsid];
  right = ts[nextid];
  news = ts[modid];
  newe = te[modid];
  
  //printf("left=%g\tright=%g\n",left,right);

  //choose move start or move end
  r = myrand(&(mempool->seed));  

  if (r < 0.5) {//shift the start
    right = te[modid];
    propose_time (&(mempool->seed),&news, left, right, beta);
  }
  else {//shift the end
    left = ts[modid];
    propose_time (&(mempool->seed),&newe, left, right, beta);
  }

  weight_shift (mempool, i_o, modid, news, newe, par);
}



void proposal_add(Mempool *mempool, Par par)
{
  double r = myrand(&(mempool->seed));
  double beta=par.beta;
  int i_o = 0;
  i_o = floor(r * 2.0 * mempool->n_orbit);
  int i;
  double news;
  double newe;
  double *ts;
  double *te;
  int nt;

  //if (mempool->isFull[i_o])
  //  return;

  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  nt = mempool->n[i_o];  


  propose_time (&(mempool->seed),&news, 0, beta, beta);  
  i=0;

  if(nt==0){
    propose_time (&(mempool->seed),&newe, 0, beta, beta);
    weight_add(mempool,i_o,i,news,newe,par);
  }
  else{

    for(i=0;i<nt;i++){

      if(ts[ (i+1) % nt]>news && te[i]<news){
	/*  ..---ts[i]~~~~~te[i]---------ts[i+1]~~~te[i+1]---  */
	/*                          ^                          */
	propose_time (&(mempool->seed),&newe, news, ts[ (i+1)% nt ], beta);
	weight_add(mempool,i_o,(i+1)%nt,news,newe,par);
	break;
      }
      if((ts[ (i+1) % nt]<te[i])&&((news<ts[ (i+1) % nt])||( news > te[i] ))){
	/*  ------ts[i+1]~~~~~te[i+1]----.....---ts[i]~~~te[i]-------  */
	/*    ^                            or                        ^     */
	propose_time (&(mempool->seed),&newe, news, ts[ (i+1) % nt], beta);
	weight_add(mempool,i_o,(i+1)%nt,news,newe,par);
	break;
      }

    }
    

  }


}





void proposal_remove(Mempool *mempool, Par par)
{

  double r = myrand(&(mempool->seed));

  int i_o = 0;

  i_o = floor(r*2*mempool->n_orbit);

  int nt = 0;

  nt = mempool->n[i_o];

  if (nt == 0)
    return;

  r = myrand(&(mempool->seed));
  int remid = floor (r * nt);

  weight_remove(mempool,i_o,remid, par);

}



void proposal_rem_anti(Mempool *mempool, Par par)
{


  double r = myrand(&(mempool->seed));

  int i_o = 0;

  i_o = floor(r*2*mempool->n_orbit);

  int nt = 0;

  nt = mempool->n[i_o];

  if (nt == 0)
    return;

  r = myrand(&(mempool->seed));
  int remid = floor (r * nt);

  weight_remove_anti(mempool,i_o,remid, par);

}


void proposal_add_anti(Mempool *mempool, Par par)
{


  double r = myrand(&(mempool->seed));

  int i_o = 0;

  i_o = floor(r*2*mempool->n_orbit);

  int nt ;
  double *ts;
  double *te;
  double beta=par.beta;
  double newe, news;

  ts = mempool->is[i_o];
  te = mempool->ie[i_o];
  nt = mempool->n[i_o];  


  //if ((mempool->isFull[i_o] ==0) && (nt == 0))
  //    return;

  r = myrand(&(mempool->seed));

  propose_time (&(mempool->seed),&news, 0, beta, beta);
  
  int i=0;


  //  if(mempool->isFull[i_o]){
  if(nt==0){
    propose_time (&(mempool->seed),&newe, 0, beta, beta);
    weight_add_anti(mempool, i_o,i,news,newe,par);
  }
  else{
    for(i=0;i<nt;i++){

      if(ts[i]<news && te[i]>news){
	propose_time (&(mempool->seed),&newe, news, te[i], beta);

	weight_add_anti(mempool, i_o ,i,news,newe,par);
	break;
      }
      if((ts[i]>te[i])&&((news>ts[i])||( news < te[i] ))){
	propose_time (&(mempool->seed),&newe, news, te[i], beta);

	weight_add_anti(mempool,i_o,i,news,newe,par);
	break;
      }

    }
    
  }


}

