#include <stdio.h>
#include <stdlib.h>



int main(){
  FILE *fp;
  
  fp=fopen("parameters.txt","r");

  int lines=0;
  int i;
  float beta,ef,u,v2;
  int iter1,iter2,iter3;
  while (EOF != (fscanf(fp,"%*[^\n]"), fscanf(fp,"%*c")))
    ++lines;

  fclose(fp);

  fp=fopen("parameters.txt","r");
  for(i=0;i<lines;i++){
    fscanf(fp,"%f\t%f\t%f\t%f\t%d\t%d\t%d\n",&beta,&ef,&u,&v2,&iter1,&iter2,&iter3);
    printf("%f\t%f\t%f\t%f\t%d\t%d\t%d\n",beta,ef,u,v2,iter1,iter2,iter3);
  }
  fclose(fp);
  printf("number of lines=%d\n",lines);
  exit(0);
  //  Par par;



}
