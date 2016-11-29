#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <assert.h>


#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "matrix.h"
#include "measurement.h"
#include "parameters.h"
#include "calc_green_func.h"
#include "utility.h"

//#include <mkl.h>

void print_average(Mempool *mem, Par par){
  //int it_average=1;
  //  printf("Preparing output...\n");

  FILE *f;

  char dir_name[100];
  sprintf(dir_name,"Job_output_%d",par.job_id);
  //  mkdir(dir_name,0777);
  char file_name[100];


  
  int task;
  int n_task=par.task;
  int n_task_i=n_task;
  int i_omega;
  int i_leg;
  int i_tau;
  int i_o;
  int n_o=mem->n_orbit*2;
  double beta=par.beta;
  int iter=0;
  double n_av=0;
  double nn_av_sum=0;
  double nn_av[N_CHANNEL] = {0};
  double err_up=0;
  double err_dn=0;
  int N_MEASURE=(par.check_err*par.n_check);

#ifdef MEASURE_GIO
  double complex green_iomega[N_ORBIT*2][N_OMEGA * 2]={0 + 0*I};
#endif

#ifdef MEASURE_LEG
  double complex green_legendre[N_ORBIT*2][N_LEG] = {0+ 0*I};
  double green_legendre2[N_ORBIT*2][N_LEG] = {0};
#endif

#ifdef MEASURE_SCPT
  double chargeScpt[N_ORBIT][N_TAU] = {0.0};
  double spinScpt[N_ORBIT][N_TAU] = {0.0};
#endif

  double complex green_tau_sum[N_ORBIT*2][N_TAU] = {0 + 0*I};


  unsigned int green_tau_total[N_ORBIT*2][N_TAU] = {0};


  int addp=0;
  int adda=0;
  int remp=0;
  int rema=0;
  int modp=0;
  int moda=0;
  int antp=0;
  int anta=0;
  int rantp=0;
  int ranta=0;
  int add_11=0;
  int add_00=0;
  int add_01=0;
  int add_10=0;
  int rem_11=0;
  /*
  for(i_o=0;i_o<2*N_ORBIT;i_o++){

#ifdef MEASURE_GIO
    for(i_omega=0;i_omega< 2 * N_OMEGA;i_omega++){
      green_iomega[i_o][i_omega] = 0.0 + 0.0 * I;
    }
#endif

#ifdef MEASURE_LEG
    for(i_leg=0;i_leg < N_LEG;i_leg++){
      green_legendre[i_o][i_leg] = 0.0 + 0.0 * I;
    }
#endif


    for (i_tau=0;i_tau<N_TAU;i_tau++){
      green_tau_sum[i_o][i_tau]=0.0;
      green_tau_total[i_o][i_tau]=0;
    }
  }
    */

  for(task=0;task<n_task;task++){
    if (mem[task].err_max < Epsilon) {
      for(i_o=0;i_o< 2 * N_ORBIT;i_o++){
	nn_av_sum+=mem[task].nn_av[i_o];
	nn_av[i_o]+=mem[task].nn_av[i_o];
	n_av+=mem[task].n_av[i_o];
      }
      addp+=mem[task].add_prop;
      adda+=mem[task].add_accpt;
      remp+=mem[task].remove_prop;
      rema+=mem[task].remove_accpt;
      modp+=mem[task].shift_prop;
      moda+=mem[task].shift_accpt;
      antp+=mem[task].add_anti_prop;
      anta+=mem[task].add_anti_accpt;
      rantp+=mem[task].rem_anti_prop;
      ranta+=mem[task].rem_anti_accpt;

      add_11+=mem[task].add_11;
      add_01+=mem[task].add_01;
      add_10+=mem[task].add_10;
      add_00+=mem[task].add_00;
      rem_11+=mem[task].rem_11;



      for(i_o=0;i_o< 2 * N_ORBIT;i_o++){
#ifdef MEASURE_GIO
#pragma simd   
	for(i_omega=0;i_omega<2*N_OMEGA;i_omega++){

#ifdef DEBUG_VERBOSE_GIO
	  printf("i_o=%d,iomega=%d,imag(sum)=%lf\n",i_o,i_omega,cimag(mem[task].green_iomega[i_o][i_omega]));
#endif
	  green_iomega[i_o][i_omega] += mem[task].green_iomega[i_o][i_omega];
	}
#endif


#ifdef MEASURE_LEG
	sprintf(file_name,"%s/green_leg%d_%d.txt",dir_name,i_o,task);  
	f=fopen(file_name,"w+");

	for(i_leg=0;i_leg< N_LEG; i_leg++){
	  double complex gl = mem[task].green_legendre[i_o][i_leg] / N_MEASURE;	  
	  double c = sqrt(2 * i_leg +1 ) /beta;
	  green_legendre[i_o][i_leg] += gl;
	  fprintf(f,"%g\n",creal(gl*c));
	  green_legendre2[i_o][i_leg] += creal(gl) * creal(gl);
	}
	fclose(f);
#endif


#pragma simd   
	for(i_tau=0;i_tau<N_TAU;i_tau++){
	  green_tau_sum[i_o][i_tau] += mem[task].green_tau[i_o][i_tau];
	  green_tau_total[i_o][i_tau]+= mem[task].green_tau_n[i_o][i_tau];
	}
      }

#ifdef MEASURE_SCPT
      for(i_o=0;i_o<N_ORBIT;i_o++)
	for(i_tau=0;i_tau< N_TAU; i_tau++){
	  chargeScpt[i_o][i_tau] += mem[task].chargeScpt[i_o][i_tau];
	  spinScpt[i_o][i_tau] += mem[task].spinScpt[i_o][i_tau];
	}
#endif

    }
    else{
      
      //      printf("%g\t%g\n",mem[task].err_up_max,mem[task].err_dn_max);
      n_task_i -- ;
    }
  }
    
  //  double correction[N_TAU*2];// Correction from tail


  int i,j;
  /*
    double tau[N_TAU*2];

#pragma simd
  for(i=0;i<2*N_TAU;i++){
    correction[i]=0;
    tau[i]=beta/N_TAU*(i+0.5-N_TAU);
  }

  for (j=0;j<N_OMEGA;j++){
    double omega=PI/beta*(2*j+1);
#pragma simd
    for(i=0;i<2*N_TAU;i++){
      correction[i]+=sin(omega*tau[i])/omega;
    }
  }

#pragma simd
  for(i=0;i<2*N_TAU;i++){
    correction[i]= correction[i]*2/beta - 0.5*((tau[i] > 0) ? 1 : -1 );
  }

  */
  //printf("Writing output...\n");

 

  for(i_o=0;i_o<2*N_ORBIT;i_o++){
#ifdef MEASURE_GIO
    sprintf(file_name,"%s/green_iomega_%d.txt",dir_name,i_o);  
    f=fopen(file_name,"w+");

    for(i_omega=0;i_omega<2*N_OMEGA;i_omega++){
      fprintf(f,"%g\t%g\n",creal(green_iomega[i_o][i_omega]/n_task_i),cimag(green_iomega[i_o][i_omega])/n_task_i);
    }

    fclose(f);
#endif

#ifdef MEASURE_LEG
    sprintf(file_name,"%s/green_legendre_%d.txt",dir_name,i_o);  
    f=fopen(file_name,"w+");
    int i_leg;
    for(i_leg=0;i_leg<N_LEG;i_leg++){
      double c = sqrt(2 * i_leg +1 ) /beta;
      //      double complex gl= - 1.0 * c * green_legendre[i_o][i_leg] / nn_av[i_o] / N_MEASURE ;
      double complex gl = c * green_legendre[i_o][i_leg] / (n_task_i);
      double gl2 = c * c * green_legendre2[i_o][i_leg];
      double error = sqrt( (gl2-creal(gl)*creal(gl) * n_task_i) / (n_task_i-1) / n_task_i);

#ifdef DEBUG_VERBOSE_LEGENDRE
      printf("c=%1.10f\n",c);
      printf("nn_av=%1.10f\n",nn_av[i_o]);
      printf("N_measure=%d\n",N_MEASURE);
      printf("%1.10f\t%1.10f\n",creal(green_legendre[i_o][i_leg]),cimag(green_legendre[i_o][i_leg]) );
      printf("gl=%1.10f+%1.10f*I\n",creal(gl),cimag(gl) );
#endif

      fprintf(f,"%1.10f\t%1.10f\t%1.10f\n",creal(gl),cimag(gl),error);
    }

    fclose(f);
#endif


    sprintf(file_name,"%s/green_tau_%d.txt",dir_name,i_o);  
    f=fopen(file_name,"w+");

    for(i_tau=0;i_tau<N_TAU;i_tau++){
      fprintf(f,"%g\t%g\t%d\n",creal(green_tau_sum[i_o][i_tau]),cimag(green_tau_sum[i_o][i_tau]), green_tau_total[i_o][i_tau]);    
    }

    fclose(f);

  }


#ifdef MEASURE_SCPT
  for(i_o=0;i_o<N_ORBIT;i_o++){
    sprintf(file_name,"%s/chargeScpt_%d.txt",dir_name,i_o);  
    f=fopen(file_name,"w+");
    for(i_tau=0;i_tau<N_TAU;i_tau++){
      fprintf(f,"%g\n",chargeScpt[i_o][i_tau]/n_task_i/N_MEASURE);
    }
    fclose(f);

    sprintf(file_name,"%s/spinScpt_%d.txt",dir_name,i_o);  
    f=fopen(file_name,"w+");
    for(i_tau=0;i_tau<N_TAU;i_tau++){
      fprintf(f,"%g\n",spinScpt[i_o][i_tau]/n_task_i/N_MEASURE);
    }
    fclose(f);
  }
#endif

  sprintf(file_name,"%s/parameters.txt",dir_name);    
  FILE *fp=fopen(file_name,"w+");

  fprintf(fp,"%f\n%f\n%f\n%f\n%f\n%d\n%d\n%d\n%d\n%d\n",par.beta,par.ef,par.u,par.v2,par.D,par.iter_measure,par.check_err,par.iter_warm,par.n_check,n_task);

  fclose(fp);


  sprintf(file_name,"%s/output.txt",dir_name);    
  fp=fopen(file_name,"w+");

  fprintf(fp,"%g\n",n_av/n_task_i);

  fprintf(fp,"%g\n%g\n%g\n%g\n%g\n",(float)adda/addp,(float)rema/remp,(float)moda/modp,(float)anta/antp,(float)ranta/rantp);
  fprintf(fp,"%g\n%g\n%g\n%g\n",(float)add_11/addp,(float)add_10/addp,(float)add_01/addp,(float)add_00/addp);

  fprintf(fp,"%d\n",n_task_i);

  fprintf(fp,"%d\n", par.proc_id);
#ifdef HAVE_MPI
  double t2= MPI_Wtime();
  fprintf(fp,"%f\n", par.t1 - par.t0);
  fprintf(fp,"%f\n", t2 - par.t0);
#endif 
  fclose(fp);

#ifdef HAVE_MPI
  sprintf(file_name,"%s/proc_name.txt",dir_name);    
  fp=fopen(file_name,"w+");


  char node_name[100];
  int resultlen;

  MPI_Get_processor_name( node_name, &resultlen );

  fprintf(fp,"processor name: %s\n", node_name);
  fclose(fp);
  //  printf("err_max=%g\n",)
#endif
  //  printf("Finish output!\n");

}


void measure_green(Mempool *mem,Par par,int iter){

  int up,i,j,sign;
  int i_o;
  double tau;
  //int task=mem->id;
  //FILE *fp1, *fp2;

  double beta=par.beta;
  int N_MEASURE=(par.check_err*par.n_check);
  
#ifdef MEASURE_GIO
  for (i = 0; i < N_OMEGA*2; i++) {
    int n = 2 * (i-N_OMEGA) + 1;
    double omega = n * PI / beta;
  
    for(i_o=0;i_o<N_ORBIT*2;i_o++){
      double complex green_s = measure_green_s(mem,i_o,beta,omega)/N_MEASURE;
      mem->green_iomega[i_o][i] += green_s;

#ifdef DEBUG_VERBOSE_GIO
      printf("i_o=%d,i=%d,omega=%g,real=%g,imag=%g\n",i_o,i,omega,creal(green_s),cimag(green_s));
      printf("i_o=%d,i=%d,sum_real=%g,sum_imag=%g\n",i_o,i,creal(mem->green_iomega[i_o][i]),cimag(mem->green_iomega[i_o][i]));
#endif
    }
    
  }
#endif

#ifdef MEASURE_LEG
  for(i_o=0;i_o<N_ORBIT*2;i_o++)
    measure_green_legendre(mem->green_legendre[i_o],mem,i_o,beta);
#endif


#ifdef MEASURE_SCPT
  for(i_o=0;i_o<N_ORBIT;i_o++)
    measure_susceptibility(mem->chargeScpt[i_o],mem->spinScpt[i_o],mem,i_o,beta);
#endif


  for(i_o=0;i_o< 2*N_ORBIT; i_o ++){
    int n=mem->n[i_o];
    for(i=0;i<n;i++){
      double e=mem->ie[i_o][i];
      for(j=0;j<n;j++){
	double s=mem->is[i_o][j];
	tau = e-s;
	sign = (tau > 0);
	tau += (1 - sign) * beta;
	int index=floor(tau/(beta/N_TAU));
	
	//      printf("%4g\t%4g\t%4g\t%8g\n",e,s,tau,creal(mem->m_up->inv[i*n+j]));
	
	mem->green_tau[i_o][index] += (2*sign-1) * mem->m[i_o]->inv[i*n+j];
	mem->green_tau_n[i_o][index]++;
      }
    }
  }

}


double complex measure_green_s(Mempool *mem,int i_o,double beta, double omega){
  int n,i,j;
  MTYPE *m;
  double *is;
  double *ie;

  n=mem->n[i_o];
  m=mem->m[i_o]->inv;
  is=mem->is[i_o];
  ie=mem->ie[i_o];
  
  double complex sum = 0.0 + 0.0 * I;
  
  
  for (i=0;i<n;i++){
    double te = ie[i];
    for(j=0;j<n;j++){
      double t=(te-is[j]);
      int sign=(t>0);
      t += (1-sign) * beta;
      t *= (omega);
      sum += ((cos(t) + I * sin(t) )* m[i*n+j]* (2*sign-1));
      
    }
  }
  sum *= (1.0/beta);
  return(sum);

  
  //size_t sz=sizeof(double complex)*n;
  //double complex *v1=(double complex*)_mm_malloc(sz,64);
  //double complex *v2=(double complex*)_mm_malloc(sz,64);

  /*  
  double complex v1[ARRAY_SZ] __attribute__((aligned(64)));
  double complex v2[ARRAY_SZ] __attribute__((aligned(64)));
  double complex v3[ARRAY_SZ] __attribute__((aligned(64)));

  int sign[ARRAY_SZ*ARRAY_SZ];

  for(i=0;i<n;i++){
    for(j=0;j<n;j++)
      sign[i*n+j]=2*(ie[i]>is[j])-1;
  }

#pragma simd
  for(i=0;i<n;i++){
    double t=ie[i]*omega;
    v1[i]=cos (t) - I * sin (t);
  }

#pragma simd
  for(i=0;i<n;i++){
    double t=is[i]*omega;
    v2[i]=cos (t) + I * sin (t);
  }
    
#pragma simd
  for(i=0;i<n;i++)
    v3[i]=0;
  
  for(j=0;j<n;j++){  
#pragma simd
    for(i=0;i<n;i++){
      v3[i]+=v2[j]*(m[i*n+j])*sign[i*n+j];
    }
  }
  
  //double complex sum=0;

  for(i=0;i<n;i++)
    sum+=v1[i]*v3[i];

  sum *= (-1.0/beta);

  */
  
  
  return sum;
  
}


void measure_green_legendre(double complex *result, Mempool *mem,int i_o,double beta){
  /*
   see arxiv:1104.3215 appendix C.2
   */

  int n,i,j,l;
  MTYPE *m;
  double *is;
  double *ie;
  MTYPE p0,p1,p2;
  MTYPE mij;
  int sign;
  double tau;
  double xt;
  MTYPE gl;
  double agl;

  n = mem->n[i_o];
  m = mem->m[i_o]->inv;
  is = mem->is[i_o];
  ie = mem->ie[i_o];

#ifdef DEBUG_VERBOSE_LEGENDRE
	printf("measurement.c:L438 n=%d\n",n);
#endif


  for (i = 0; i < n; i ++){
    for(j = 0; j < n; j ++){
      mij = m[i * n + j];
      tau = ie[i] - is[j];
      sign = 2 * (tau >= 0) -1;
      tau += (tau < 0) * beta;
      xt = 2.0 * tau / beta -1.0;
      gl = mij*sign;
      result[0] += gl;

#ifdef DEBUG_VERBOSE_LEGENDRE
      printf("%1.10f\t%1.10f\n",creal(result[0]),cimag(result[0]) );
#endif      

      gl = gl*xt;
      result[1] += gl;
      p0 = 1.0;
      p1 = xt;
      for(l = 2; l < N_LEG; l ++){
	p2 = ((2.0 * l - 1.0) * xt * p1 - (l - 1.0) * p0) / l;
	gl = mij * p2 * sign;
	result[l] += gl;	
	
	p0 = p1;
	p1 = p2;
      }
    }
  }
}


void measure_susceptibility(double *chargeScpt, double * spinScpt, Mempool *mem,int i_o, double beta){
  double dt;
  int iSeg,iTau;
  double tStart,tEnd;

  /* n_up(tau)n_up(0),  n_dn(tau)n_dn(0),  n_up(tau)n_dn(0),  n_dn(tau)n_up(0) */
  double ouu,odd,oud,odu;

  int nu = mem->n[i_o*2];
  int nd = mem->n[i_o*2+1];
  double * isu = mem->is[i_o*2];
  double * ieu = mem->ie[i_o*2];
  double * isd = mem->is[i_o*2+1];
  double * ied = mem->ie[i_o*2+1];

  for(iTau=0;iTau<N_TAU;iTau++){
    dt = iTau*beta/N_TAU;
    ouu = 0.0;
    oud = 0.0;
    odu = 0.0;
    odd = 0.0;

    for(iSeg=0;iSeg<nu;iSeg++){
      tStart=isu[iSeg]+dt;
      tEnd=ieu[iSeg]+dt;
      tStart -= (tStart > beta) * beta;
      tEnd -= (tEnd > beta) * beta;
      ouu += overlap_total(isu,ieu, nu, tStart, tEnd, beta);
      oud += overlap_total(isd,ied, nd, tStart, tEnd, beta);   
    }
    for(iSeg=0;iSeg<nd;iSeg++){
      tStart=isd[iSeg]+dt;
      tEnd=ied[iSeg]+dt;
      tStart -= (tStart > beta) * beta;
      tEnd -= (tEnd > beta) * beta;
      odu += overlap_total(isu,ieu, nu, tStart, tEnd, beta);
      odd += overlap_total(isd,ied, nd, tStart, tEnd, beta);   
    }
    chargeScpt[iTau] += (ouu+oud+odu+odd);
    spinScpt[iTau] += (ouu-oud-odu+odd);
  }
}
