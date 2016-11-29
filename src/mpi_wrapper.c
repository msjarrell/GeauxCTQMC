#ifdef HAVE_MPI
// mpicc -std=c99  a.c

// alternative implemenation:
//    atomic imcrement on global variable


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <errno.h>
#include <mpi.h>
#include <unistd.h>


#include <sys/stat.h>
#include <sys/types.h>

#include "matrix.h"
#include "main.h"

#define QUEUE_SIZE 100
#define MANAGER_ID 0
#define FINISH_SIGNAL INT32_MAX

typedef struct
{
  int id;
} Job;

/*
void Myleep (int us)
{
    struct timespec t;
    t.tv_sec = us / (1000 * 1000);
    t.tv_nsec = (us % (1000 * 1000)) * 1000;
    nanosleep (&t, NULL);
}
*/



void InitQueue (Job *queue,  int queue_size)
{
  //srand ((long) time (NULL))
  //Par *pars=(Par *)malloc(sizeof(Par)*queue_size);
  //  init_CTQMC(pars,queue_size);

  int i;
  for (i=0;i<queue_size;i++)
    queue[i].id=i;
}


void DoSomething (const Job job, const int id,double t0)
{
  //MySleep ((rand () % 500 * 1000));
  //sleep ((rand () % 3));
  printf ("worker %02d starting job %05d\n", id, job.id);
  /*
  sleep (job.parameter);
  char dir_name[100];
  sprintf(dir_name,"Job_output_%d",job.id);
  mkdir(dir_name,0777);
  */


  main_CTQMC(job.id, t0, id);

  printf ("worker %02d finished job %05d\n", id, job.id);
}


void Manager (const int id, Job *queue, int queue_size, const int nworker)
{

  //  Job *queue = (Job *) malloc(sizeof (Job) * queue_size);
  InitQueue (queue, queue_size);
  Job *queue_ptr = queue; // queue pointer, where a new job is dispatched

  int i;
  int recv_msg[1];
  int send_msg[1];
  MPI_Status status;
  const int mytag = 0;


  printf("Manager starts working\n");
  while (queue_ptr < queue + queue_size) {
    MPI_Recv (recv_msg, 1, MPI_INT, MPI_ANY_SOURCE, mytag, MPI_COMM_WORLD, &status);
    send_msg[0] = queue_ptr->id;
    //    send_msg[1] = queue_ptr->par;
    queue_ptr++;
    const int dest = recv_msg[0];
    MPI_Send (send_msg, 1, MPI_INT, dest, mytag, MPI_COMM_WORLD);
  }


  // send finish signals
  for (i = 0; i < nworker; ++i) {
    MPI_Recv (recv_msg, 1, MPI_INT, MPI_ANY_SOURCE, mytag, MPI_COMM_WORLD, &status);
    send_msg[0] = FINISH_SIGNAL;
    //send_msg[1] = FINISH_SIGNAL;
    const int dest = recv_msg[0];
    MPI_Send (send_msg, 1, MPI_INT, dest, mytag, MPI_COMM_WORLD);
    printf ("\t\t\t\t manager: kill worker %02d\n", dest);
  }



  free (queue);

}



void Worker (const int id)
{
  int send_msg[1] = {id};
  int recv_msg[2];
  const int dest = MANAGER_ID;
  MPI_Status status;
  const int mytag = 0;

  printf("Worker %d started working\n", id);

  double t0;
  t0 = MPI_Wtime();


  while (1) {
    MPI_Send (send_msg, 1, MPI_INT, dest, mytag, MPI_COMM_WORLD);
    MPI_Recv (recv_msg, 1, MPI_INT, MPI_ANY_SOURCE, mytag, MPI_COMM_WORLD, &status);

    if (recv_msg[0] == FINISH_SIGNAL ) {
      printf ("\t\t\t\t\t\t\t worker %02d retired\n", id);
      break;
    }
    else {
      Job job;
      job.id = recv_msg[0];
      //      job.parameter = recv_msg[1];
      DoSomething (job, id, t0);
    }
  }

}




int main (int argc, char *argv[])
{
  int id, nprocs;
  Job* queue;
  



  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);



  if(nprocs<2){
    printf("Need at least one worker to start job!\n");
    exit(0);

  }
  /*
  char name[100];
  int resultlen;
  MPI_Get_processor_name( name, &resultlen );
  printf("%d: %s\n",id, name);
  */
    



  int queue_size=get_sz();
  
  queue = (Job *) malloc(sizeof (Job) * queue_size);
  
  if (id == MANAGER_ID)
    Manager (id, queue,queue_size, nprocs - 1);
  else
    Worker (id);

  MPI_Finalize();


  return 0;
}

#endif
