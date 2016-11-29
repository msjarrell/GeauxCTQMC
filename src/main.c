#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>

#include "complex.h"

#include "matrix.h"
#include "parameters.h"
#include "main.h"

#ifndef USE_MATRIX
#include "weight.h"
#endif

#ifdef USE_PCG
#include "pcg-c-basic-0.9/pcg_basic.h"
#endif

#ifdef USE_MATRIX
#include "weight_dhm.h"
#endif

#include "utility.h"
#include "calc_green_func.h"
#include "measurement.h"
#include "proposal.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <omp.h>


int get_sz(){
  FILE *fp;
  fp=fopen("parameters.txt","r");

  int lines=0;
  int i;
  while (EOF != (fscanf(fp,"%*[^\n]"), fscanf(fp,"%*c")))
    ++lines;
  fclose(fp);
  return(lines);
}


#ifndef HAVE_MPI
void main(){
  int jobs=get_sz();
  //Par *pars=(Par *)malloc(sizeof(Par)*jobs);
  //  init_CTQMC(pars,jobs);
  int i;
  for(i=0;i<jobs;i++){
#ifdef VERBOSE
    printf("Executing job %d\n",i);
#endif
    main_CTQMC(i,0,0);
#ifdef VERBOSE
    printf("Finishing job %d\n",i);
#endif
    }

  //  free(pars);
  exit(0);
}
#endif


int init_CTQMC(Par *par,int job_id,double t0,int proc_id){
  //read in parameters from line %job_id in parameters.txt
  

  FILE *fp;
  //int jobs=queue_size;
  //  Par *pars=(Par *)malloc(sizeof(Par)*lines);

  float beta,ef,u,v2,D;
  int iter1,iter2,iter3,iter4;
  int i;
  int task;
  int seed;
  //printf("Proc %d Loading parameters for job %d\n",proc_id,job_id);

  fp=fopen("parameters.txt","r");

  i=0;
  while(i<job_id){
    fscanf(fp,"%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n",&beta,&ef,&u,&v2,&D,&iter1,&iter2,&iter3,&iter4,&task,&seed);
    i++;
  }

  fscanf(fp,"%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n",&beta,&ef,&u,&v2,&D,&iter1,&iter2,&iter3,&iter4,&task,&seed);
  par->beta=beta;
  par->u = u;
  par->ef = ef;
  par->v2 = v2;
  par->D = D;
  par->iter_measure = iter1;
  par->check_err = iter2;
  par->n_check  = iter3;
  par->iter_warm  = iter4;
  par->task=task;
  par->job_id=i;
  par->proc_id=proc_id;
  par->seed=seed;

  fclose(fp);

  par->t0=t0; 

#ifdef HAVE_MPI
  par->t1= MPI_Wtime();
#endif
#ifdef VERBOSE
  printf("Load parameters finished\n");
#endif
  char dir_name[100];
  sprintf(dir_name,"Job_output_%d",par->job_id);
  mkdir(dir_name,0777);

  return(0);
  //exit(0);
}


int
main_CTQMC (int job_id, double t0, int proc_id)
{
  Par par;
  
 
  init_CTQMC(&par, job_id,t0, proc_id);

#ifdef USE_MATRIX
  //read the hamiltonian
  //  printf("load hamiltonian\n");
  init_hamiltonian(&par);
  //  printf("load hamiltonian finished\n");
#endif
  
  init_g_tau (par.g_tau,par.beta,par.v2, par.D);	//initialize the table for green function
  

  if(0)
    {
      //Use single thread
      return main_serial(par);
    }
  else
    {
      //use OpenMP
      return main_parallel(par);
    }
  
  return(0);
}

int
main_serial ( Par parameters)
{
  //msrand(time(NULL));
  //seed[id]= time(NULL) ^ id;

  Mempool mempool;

  mem_alloc(&mempool);

#ifdef DEBUG_VERBOSE
  printf("Mem allocated\n");
#endif

  if (parameters.seed == -1)
#ifdef USE_PCG
    pcg32_srandom_r(&(mempool.seed), time(NULL) ^ (intptr_t)&printf, 0);
#else
    mempool.seed=time(NULL);
#endif
    else
#ifdef USE_PCG
    pcg32_srandom_r(&(mempool.seed), parameters.seed, 0);
#else
    mempool.seed = parameters.seed;
#endif
  mempool.id=0;
  main_loop (&mempool,parameters);
  mem_free(&mempool);
  return 0;
}


int
main_parallel (Par par)
{

  //  Mempool mempool[TASK];
  int task=par.task;
  Mempool *mempool=(Mempool *)malloc(sizeof(Mempool)*task);
  if (mempool == NULL){
    printf("memory allocation failed!\n");
    exit(0);
  }
  int i;

#ifdef DEBUG_VERBOSE
  printf("Trying weight_loc\n");
#endif


#ifdef USE_MATRIX
  double trace=weight_loc(0,NULL,NULL, par);
#endif

  for(i=0;i<task;i++){

    mempool[i].n_orbit=N_ORBIT;

#ifdef DEBUG_VERBOSE
    printf("Mem getting allocated i=%d\n",i);
#endif

    mem_alloc(mempool+i);
    mempool[i].id=i;

#ifdef DEBUG_VERBOSE
    printf("Mem allocated\n");
#endif


#ifdef USE_MATRIX
    mempool[i].trace=trace;
#endif

#ifdef USE_PCG
    if (par.seed == -1)
      pcg32_srandom_r(&(mempool[i].seed), time(NULL) ^ (intptr_t)&printf, i);
    else
      pcg32_srandom_r(&(mempool[i].seed), par.seed, i);

#else
    if (par.seed == -1)
      mempool[i].seed = time(NULL) ^ i;
    else
      mempool[i].seed = par.seed ^ i;

#endif
  }

#pragma omp parallel for 
  for(i=0;i<task;i++)   {  
      main_loop(mempool+i, par);
      //printf("TASK %3d/%3d Finished\n",i,TASK);
  }


  print_average(mempool, par);
  //  measure_green_f(mempool,par);

  for(i=0;i<task;i++)
    mem_free(mempool+i);


  free(mempool);
    //printf("----------------------\n");
  //printf("Total time=%f\n",total_time);
  //MPI_Finalize ();
  return 0;
}


void
main_loop (Mempool *mempool,Par par)
{

  //main_init (mempool, g_tau, BETA);	//put one segment on each channel to start

  //print_seg(ius,iue,nu,ids,ide,nd);

  int i,j,k;


  int N_MEASURE=(par.check_err*par.n_check);
  int n_o=mempool->n_orbit*2;
  
  double news,newe;
  double beta=par.beta;


  for (i = 0; i < par.iter_warm; i++) {
    for(j=0;j< par.check_err;j++){
      for(k=0;k< par.iter_measure;k++){
	proposal_next (mempool, par);
      }
    }

    for (j=0;j<n_o;j++){
      double error;
      error=inverse(mempool->m[j]);
      if(error > mempool->err_max)
	mempool->err_max = error;	
    }
  }
  for (i=0; i < par.n_check; i++) {
    for(j=0;j < par.check_err;j++){
      for(k=0;k < par.iter_measure;k++){
	proposal_next (mempool, par);
      }
      measure_green(mempool,par, i*par.check_err+j);
      for(k=0;k<n_o;k++){
	
	double nOp = (double)mempool->n[k];
#ifdef DEBUG_VERBOSE_LEGENDRE
	printf("%d\t%lf\n",mempool->n[k],nOp);
#endif
	mempool->n_av[k] += nOp / N_MEASURE;  
	mempool->nn_av[k] += nOp * nOp / N_MEASURE;  
#ifdef DEBUG_VERBOSE_LEGENDRE
	printf("n_av=%lf\n",mempool->n_av[k]);
	printf("nn_av=%lf\n",mempool->nn_av[k]);
#endif
      }
    }
    for(j=0;j<n_o;j++){
      double error=inverse(mempool->m[j]);    
      if(error > mempool->err_max)
	mempool->err_max = error;	
    }

  }

}


void __attribute__((optimize("O0")))
mem_alloc (Mempool *mem){
  
  int n_o= mem->n_orbit*2;
  int i,j;
  //  printf("allocating mem\n");


  mem->m_tmp1=(Matrix *)malloc(sizeof(Matrix));
  mem->m_tmp2=(Matrix *)malloc(sizeof(Matrix));

#ifdef DEBUG_VERBOSE
  printf("starting allocating n_o=%d\n",n_o);
#endif
  size_t msize= sizeof(Matrix *) * n_o;
  //  printf("msize=%d\n",msize);
  void *ptr = malloc(msize);
#ifdef DEBUG_VERBOSE
  printf("b1\n");
  printf("test\n");
  printf("allocating mem->m\n");
  
#endif

  for(i=0;i<n_o;i++)
    mem->m[i]=(Matrix *)malloc(sizeof(Matrix));

#ifdef DEBUG_VERBOSE
  printf("allocating mem->m[i]\n");
#endif

  size_t matrix_sz=sizeof(MTYPE) *ARRAY_SZ *ARRAY_SZ;
  size_t array_sz=sizeof(double) *ARRAY_SZ;
  size_t mkl_sz=sizeof(MKLTYPE) *ARRAY_SZ *ARRAY_SZ;
  
  for(i=0;i<n_o;i++){
    mem->m[i]->mat=(MTYPE *)_mm_malloc(matrix_sz,64);
    mem->m[i]->inv=(MTYPE *)_mm_malloc(matrix_sz,64);
    mem->m[i]->mkl=(MKLTYPE *)malloc(mkl_sz);
    mem->m[i]->N=0;
  }

#ifdef DEBUG_VERBOSE
  printf("allocating mem->m_tmp\n");
#endif


  mem->m_tmp1->mat=(MTYPE *)_mm_malloc(matrix_sz,64);
  mem->m_tmp1->inv=(MTYPE *)_mm_malloc(matrix_sz,64);
  mem->m_tmp1->mkl=(MKLTYPE *)malloc(mkl_sz);
  mem->m_tmp1->N=0;

  mem->m_tmp2->mat=(MTYPE *)_mm_malloc(matrix_sz,64);
  mem->m_tmp2->inv=(MTYPE *)_mm_malloc(matrix_sz,64);
  mem->m_tmp2->mkl=(MKLTYPE *)malloc(mkl_sz);
  mem->m_tmp2->N=0;


  //mem->is=(double **)malloc(sizeof(double *)*n_o);
  //mem->ie=(double **)malloc(sizeof(double *)*n_o);

#ifdef DEBUG_VERBOSE
  printf("initializing n\n");
#endif

  for(i=0;i < n_o;i++){
    //mem->is[i]=(double *)_mm_malloc(array_sz,64);
    //mem->ie[i]=(double *)_mm_malloc(array_sz,64);
#ifdef DEBUG_VERBOSE
    printf("i=%d\n",i);
#endif

    mem->n[i] = 0;
    mem->n_av[i] = 0.0;
    mem->nn_av[i] = 0.0;
    mem->isFull[i] = 0;
   }

#ifdef DEBUG_VERBOSE
  printf("initializing n_sum\n");
#endif

  mem->n_sum=0;

#ifdef DEBUG_VERBOSE
  printf("allocating op_lists\n");
#endif



#ifdef USE_MATRIX
  mem->op_list=(int *)malloc(sizeof(int)*ARRAY_SZ*4);
  mem->tmp_op_list=(int *)malloc(sizeof(int)*ARRAY_SZ*4);

  mem->time_list=(double *)malloc(sizeof(double)*ARRAY_SZ*4);
  mem->tmp_time_list=(double *)malloc(sizeof(double)*ARRAY_SZ*4);
#endif



  mem->err_max=0.0;
  mem->add_prop=0;
  mem->add_accpt=0;
  mem->remove_prop=0;
  mem->remove_accpt=0;
  mem->shift_prop=0;
  mem->shift_accpt=0;
  mem->add_anti_prop=0;
  mem->add_anti_accpt=0;
  mem->rem_anti_prop=0;
  mem->rem_anti_accpt=0;

  mem->add_11=0;
  mem->add_10=0;
  mem->add_01=0;
  mem->add_00=0;

  mem->rem_11=0;
  mem->rem_10=0;
  mem->rem_01=0;
  mem->rem_00=0;


  mem->iter=0;

#ifdef DEBUG_VERBOSE
  printf("allocating measurements\n");
#endif

  for(i=0;i< N_CHANNEL;i++){

#ifdef DEBUG_VERBOSE
  printf("allocating g_iw\n");
#endif

#ifdef MEASURE_GIO
    for (j = 0; j < 2 * N_OMEGA; j++) {
      mem->green_iomega[i][j] = 0.0 + 0.0* I;
    }
#endif

#ifdef DEBUG_VERBOSE
  printf("allocating g_leg\n");
#endif

#ifdef MEASURE_LEG
    for (j = 0; j < N_LEG; j++) {
      mem->green_legendre[i][j] = 0.0 + 0.0 * I;
    }
#endif


#ifdef DEBUG_VERBOSE
    printf("allocating g_tau\n");
#endif

    for(j=0;j<N_TAU;j++){
#ifdef DEBUG_VERBOSE
      printf("j=%d\n",j);
#endif
      mem->green_tau[i][j] = 0.0 + 0.0 * I;
      mem->green_tau_n[i][j] = 0;
    }
#ifdef DEBUG_VERBOSE
    printf("allocating g_tau DONE\n");
#endif


  }

#ifdef MEASURE_SCPT
  for (i=0; i < N_ORBIT; i++)
    for (j = 0; j < N_TAU; j++) {
      mem->chargeScpt[i][j]= 0.0;
      mem->spinScpt[i][j]= 0.0;
    }
#endif



#ifdef DEBUG_VERBOSE
  printf("allocating mem finished\n");
  printf("one more line...\n");
#endif
}


void
mem_free (Mempool *mem){
  int i;
  int n_o= mem->n_orbit*2;

  for(i=0;i<n_o;i++){
    _mm_free(mem->m[i]->mat); 
    _mm_free(mem->m[i]->inv);
    free(mem->m[i]->mkl);
    free(mem->m[i]);
  }

  _mm_free(mem->m_tmp1->mat);
  _mm_free(mem->m_tmp1->inv);
  free(mem->m_tmp1->mkl);

  _mm_free(mem->m_tmp2->mat);
  _mm_free(mem->m_tmp2->inv);
  free(mem->m_tmp2->mkl);

  //  free(mem->m);
  free(mem->m_tmp1);
  free(mem->m_tmp2);

#ifdef USE_MATRIX
  free(mem->op_list);
  free(mem->tmp_op_list);
  
  free(mem->time_list);
  free(mem->tmp_time_list);
#endif

}

