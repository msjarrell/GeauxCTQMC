#include <cstdio>
#include <iostream>
#include <new>

void DIA2SDF(double **DIA_D, int *DIA_Did, int N, int **SDF)
// this version of DIA2SDF ignores the values
// SDF only returns the structure on NNZs using Stridden Diagon Format
// SDF should have length NNZD on allocation before being passed into here
{

  int offset;
  int base;
  int depth;
  int SV[8]={0,0,0,0,0,0,0,0};//the array length of 8 is arbitrary
  int SL[8]={0,0,0,0,0,0,0,0};//for N=1024, 5/4 is enough for SV/SL
  int counter[8]={0,0,0,0,0,0,0,0};
  int zeroes[8][N];
  int values[N];
  int k=0;
  int m=0;
  int NNZD=0;
  int bound=0;

  // this code can be used to calculate NNZD to allocate **SDF outside
  // int NNZD=0;
  // for(int i=0; i<N; i++){
  //   if(DIA_Did[i]!=0.0){
  //     NNZD++;
  //   }
  // }


  for(int i=0; i<N; i++){
    for (int ii=0; ii<8; ii++){
      SV[ii]=0;
      SL[ii]=0;
      for(int j=0; j<N; j++){
	zeroes[ii][j]=0;
      }
    }
    if(DIA_Did[i]!=0){
      NNZD++;
      for (int ii=0; ii<8; ii++){
	counter[ii]=0;
      }
      offset=i;
      k=0;
      for(int j=0; j<i+1; j++){
	if(DIA_D[i][j]==0.0){
	  counter[0]++;
	} else {
	  zeroes[0][k]=counter[0];
	  counter[0]=0;
	  k++;
	}
      }
      zeroes[0][k]=counter[0];
      counter[0]=0;
      k++;
      depth=1;
      base=zeroes[0][0];
      values[0]=zeroes[0][1];
      for(int j=2; j<k-1; j++){
	for(int jj=0; jj<depth; jj++){
	  if(values[jj]==zeroes[0][j]){
	    goto next_value;
	  }
	}
	values[depth]=zeroes[0][j];
	depth++;
      next_value:
	;
	// Again, I have assumed that when starting a new stridden
	// diagonal, after the initial run of zeroes (base), we'll start
	// from the smallest stride value and that we won't start from the middle
	// in the generic SDF, this does not necessarily hold!
	// FIX: make a historgram of each different value in zeroes and sort values accordingly
      }
      // this section is to check and see if the last stride of zeroes needs to be accounted for separately
      m=0;
      for(int jj=0; jj<depth; jj++){
	if(values[jj]>zeroes[0][k-1]){
	  m=1;
	}
      }
      if(!m){
	values[depth]=zeroes[0][k-1];
	depth++;
      }
      SV[0]=values[0]+1;
      if(depth>1){
	while(zeroes[0][SL[0]+1]==values[0]){
	  SL[0]++;
	}
	SL[0]++;
	m=0;
	for(int jj=1; jj<k; jj++){
	  if(zeroes[0][jj]!=values[0]){
	    zeroes[1][m]=zeroes[0][jj];
	    m++;
	  }
	}
	bound=SV[0]*(SL[0]-1);
	SV[1]=values[1]+1+bound;
	if(depth>2){
	  for(int j=1; j<depth-1; j++){
	    while(zeroes[j][SL[j]]==values[j]){
	      SL[j]++;
	    }
	    SL[j]++;
	    k=m;
	    m=0;
	    for(int jj=0; jj<k; jj++){
	      if(zeroes[j][jj]!=values[j]){
		zeroes[j+1][m]=zeroes[j][jj];
		m++;
	      }
	    }
	    bound+=SV[j]*(SL[j]-1);
	    SV[j+1]=values[j+1]+1+bound;
	  }
	}
      }
      SDF[NNZD-1]=new int[2*depth+2];
      SDF[NNZD-1][0]=offset-(N-1);
      SDF[NNZD-1][1]=base;
      SDF[NNZD-1][2]=depth;
      for(int l=0; l<depth; l++){
  	SDF[NNZD-1][l+3]=SV[depth-l-1];
      }
      for(int l=0; l<depth-1; l++){
  	SDF[NNZD-1][3+depth+l]=SL[depth-l-2];
      }
    }
  }
  return;
}
