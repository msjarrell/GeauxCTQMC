#ifdef HAVE_MPI
// mpicc -std=c99  a.c

// alternative implemenation:
//    atomic imcrement on global variable


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <errno.h>
//#include <mpi.h>
#include <unistd.h>


#include <sys/stat.h>
#include <sys/types.h>

#include "matrix.h"
#include "main.h"





typedef struct
{
  int id;
} Job;


int main_CTQMC (int job_id, double t0, int proc_id);
int get_sz();

void InitQueue (Job *queue,  int queue_size)
{
  //srand ((long) time (NULL))
  //Par *pars=(Par *)malloc(sizeof(Par)*queue_size);
  //  init_CTQMC(pars,queue_size);

  int i;
  for (i=0;i<queue_size;i++)
    queue[i].id=i;
}


void DoSomething (const Job job, const int id)
{
  printf ("worker %02d starting job %05d\n", id, job.id);
  main_CTQMC(job.id,0,0);
  printf ("worker %02d finished job %05d\n", id, job.id);
}



int main (int argc, char *argv[])
{
  int id, nprocs;
  Job* queue;
  int i;
    
  int queue_size=get_sz();
  
  queue = (Job *) malloc(sizeof (Job) * queue_size);

  InitQueue (queue, queue_size);
  
  for (i=0;i<queue_size;i++)
    DoSomething (queue[i], 0);
  
  free(queue);
  return 0;
}

#endif
