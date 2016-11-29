#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
  double e0;
  double g;
  double w0;
  double a;
  double U0;
} Par_H;

typedef struct{
  int n;
  double *m;
} Matrix_H;

Par_H input(){
  FILE *fp;
  double e0;
  double g;
  double w0;
  double a;
  double U0;
  Par_H par;
  
  fp=fopen("parameters_H.txt","r");
  fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf",&e0,&g,&w0,&a,&U0);
  fclose(fp);
  par=(Par_H){e0,g,w0,a,U0};
  return par;
}

Matrix_H init_matrix_h(){
  Matrix_H mh;
  int i;
  int nn=8*8;
  mh.n=8;
  mh.m=(double *)malloc(sizeof(double)*nn);

  for (i=0;i<nn;i++)
    mh.m[i]=0.0;

  return mh;
}




int main(){
  int i;
  int j;
  double e0;
  double g;
  double w0;
  double a;
  double U0;

  Par_H par= input();
  Matrix_H mh=init_matrix_h();

  e0 = par. e0 ; 
  g  = par. g  ;
  w0 = par. w0 ;
  a  = par. a  ;
  U0 = par. U0 ;

  double ee= w0 * sqrt(1+g*g);
  mh.m[0] = ee;
  mh.m[9] = ee + e0 + g*w0;
  mh.m[18] = ee + e0 + g*w0;
  mh.m[27] = ee + 2.0*e0 + U0 - (a-1.0)*g*w0;

  mh.m[4] = ee;
  mh.m[13] = ee + e0 + w0;
  mh.m[22] = ee + e0 + w0;
  mh.m[31] = ee + 2*e0 + U0;

  mh.m[32] = ee;
  mh.m[41] = ee + e0 + w0;
  mh.m[50] = ee + e0 + w0;
  mh.m[59] = ee + 2*e0 + U0;

  mh.m[36] = ee;
  mh.m[45] = ee + e0 - g*w0;
  mh.m[54] = ee + e0 - g*w0;
  mh.m[63] = ee + 2.0*e0 + U0 + (a-1.0)*g*w0;


  /*
  for(i=0;i<8;i++){
    for(j=0;j<8;j++)
      printf("\t%lf",mh.m[i*8+j]);
    printf("\n");
  }
  */
  
  for(i=0;i<8;i++){
    for(j=0;j<8;j++)
      if (mh.m[i*8+j]!=0)
	printf("%lf\t%d\t%d\n",mh.m[i*8+j],i,j);
    //    printf("\n");
  }

  return 0;
}
