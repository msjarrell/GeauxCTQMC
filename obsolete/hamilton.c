#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>

#define SITES 7
#define STATES 16384//=4^SITES

#define LOOPS 1

#define U 1.0
#define J 0.1
#define Jx (J)
#define Jp (J)
#define Up (U-2.0*J)

typedef struct{
  int n;//one or zero
  unsigned int state;
}Ket;


void diag(double *H){

  MKL_INT info;

  //int m=0;
  double w[STATES];
  //double z[STATES*STATES];
  //double *z=malloc(sizeof(double)*STATES*STATES);
  //int isupz[STATES];
  //info=LAPACKE_dsyevr(LAPACK_ROW_MAJOR,'V','A','U',STATES,H,STATES,0,1,1,10,-1.0,&m,w,z,STATES,isupz);

{ 
  info=LAPACKE_dsyevd(LAPACK_ROW_MAJOR,'V','U',STATES,H,STATES,w);
}

  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }

  //int i,j;
  //printf("m=%d\n",m);
  /*
  for(i=0;i<STATES;i++){
    printf("%d,%.3f:\t",i,w[i]);
    for(j=0;j<STATES;j++){
      printf("%.3f\t",H[j*STATES+i]);
    }
    printf("\n");

  }
  */

  //free(z);
}


void printket(Ket phi){
  int i;
  printf("%d x ",phi.n);
  for(i=0;i<SITES*2;i++){
    printf("%d",(phi.state >> i) & 0x1);
  }
  printf("\n");
}


Ket Nms(int m, int s, Ket phi){
  Ket new=phi;
  int shift=m*2+s;
  new.n=phi.n*((phi.state>>shift)&0x1);
  return new;
}

Ket dms(int m, int s, Ket phi){
  Ket new=phi;
  int shift=m*2+s;
  new.n=phi.n*((phi.state>>shift)&0x1);
  new.state=phi.state^(0x1<<shift);
  return new;
}

Ket ddms(int m, int s, Ket phi){
  Ket new=phi;
  int shift=m*2+s;
  new.n=phi.n*(1-(phi.state>>shift)&0x1);
  new.state=phi.state^(0x1<<shift);
  return new;
}




int Eq(Ket bra,Ket ket){
  int n=1;
  n=(bra.state==ket.state)*ket.n;
  return n;
}

void main(){

  
  int i,j,m,n,s;
 
  Ket ket;
  ket.state=5;
  ket.n=1;
  
  Ket bra;

  //double H[STATES*STATES];
  //double tmp[STATES*STATES];
  double *H=malloc(sizeof(double)*STATES*STATES);
  double *tmp=malloc(sizeof(double)*STATES*STATES*LOOPS);
  
  for(i=0;i<STATES;i++){
    for(j=0;j<STATES;j++){
      
      bra.state=j;
      ket.state=i;
      
      ket.n=1;
      
      
      double h=0.0;
      //      if(i<=j){
      for(m=0;m<SITES;m++){

	h+=Eq(bra,Nms(m,1,Nms(m,0,ket)))*U;
	    
	for(n=0;n<SITES;n++){
	  if(m>n){
	    for(s=0;s<2;s++){
	      h+=Eq(bra,Nms(m,s,Nms(n,1-s,ket)))*Up+Eq(bra,Nms(m,s,Nms(n,s,ket)))*(Up-J);
	    }
	  }

	  if(m>n){
	    h-=Eq(bra,ddms(m,1,ddms(n,0,dms(n,1,dms(m,0,ket)))))*Jx;
	    h-=Eq(bra,ddms(m,0,ddms(n,1,dms(n,0,dms(m,1,ket)))))*Jx;
	    h+=Eq(bra,ddms(n,0,ddms(n,1,dms(m,0,dms(m,1,ket)))))*Jp;
	    h+=Eq(bra,ddms(m,0,ddms(m,1,dms(n,0,dms(n,1,ket)))))*Jp;
	  }
	  
	  
	}
      }
      //      }      
      H[i*STATES+j]=h;

    }




  }
  
  
    
  for(i=0;i<STATES;i++){
    for(j=0;j<STATES;j++)
      if(H[i*STATES+j]!=0)
	printf("%.3f\t%d\t%d\n",H[i*STATES+j],i,j);
      //printf("%.3f\t",H[i*STATES+j]);
      
      //printf("\n");
  }
  
  /*
  for(i=0;i<LOOPS;i++){
    memcpy(tmp+i*STATES*STATES,H,sizeof(double)*STATES*STATES);
  }

#pragma omp parallel for
  for(i=0;i<LOOPS;i++)
    diag(tmp+i*STATES*STATES);

  free(H);
  free(tmp);
*/
}
