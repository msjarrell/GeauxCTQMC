#include <stdlib.h>
#include <stdio.h>

//fwds

/*void accept_modify (double **ius, double **iue, double **ids, double **ide,
		    int *nu, int *nd, int up, int modid, double news,
		    double newe);


*/
//implementations

void
accept_shift (double *ts, double *te,  int modid, double news, double newe)
{
  ts[modid] = news;
  te[modid] = newe;
}


void accept_add(double *ts,double *te,int *nt, int newid,double news,double newe){
  int i;
  for (i=*nt-1;i>=newid;i--){
    ts[i+1]=ts[i];
    te[i+1]=te[i];
  }
  ts[newid]=news;
  te[newid]=newe;
  *nt=*nt+1;
}


void accept_add_anti(double *ts,double *te,int *nt, int newid,double news,double newe){
  int i;
  if (*nt ==0){
    *nt=1;
    ts[0]=newe;
    te[0]=news;
    return;
  }
  for (i=*nt-1;i>=newid;i--){
    ts[i+1]=ts[i];
    te[i+1]=te[i];
  }
  ts[newid+1]=newe;
  te[newid]=news;
  *nt=*nt+1;
}

void accept_remove(double *ts,double *te,int *nt, int remid){
  int i;
  for (i=remid;i<*nt-1;i++){
    ts[i]=ts[i+1];
    te[i]=te[i+1];
  }
  *nt=*nt-1;
}

void accept_remove_anti(double *ts,double *te,int *nt, int remid){
  int i;
  int k=*nt;
  if (k==1){
    *nt=0;
    return;
  }

  if(remid==(k-1)){
    te[k-1]=te[0];
    for(i=0;i<k-1;i++){
      ts[i]=ts[i+1];
      te[i]=te[i+1];
    }
  }
  else{
    te[remid]=te[remid+1];
    for (i=remid+1;i<k-1;i++){
      ts[i]=ts[i+1];
      te[i]=te[i+1];
    }
  }
  *nt= k-1;
}
