double inverse(Matrix *m);

void print_seg (Mempool *mempool);
void print_matrix(Matrix *m);



void cvt_to_MKL(MTYPE *a,MKLTYPE *m, int n);

void cvt_from_MKL(MTYPE *a,MKLTYPE *m, int n);
/*
void det_update_fast(Matrix *m, Matrix *temp, MTYPE *u,MTYPE *v);

void det_add_fast (MTYPE *storage, double *ts, double *te, int nt, int newid,
		   double news, double newe, Matrix * m, Matrix * temp1, Matrix *temp2,
	      Par par);

void det_add_anti_fast (MTYPE *storage, double *ts, double *te, int nt, int newid,
		   double news, double newe, Matrix * m, Matrix * temp1, Matrix *temp2,
	      Par par);


void det_remove_fast (MTYPE *storage, double *ts, double *te, int nt, int remid,
		      Matrix * m, Matrix * temp1, Matrix *temp2, Par par);

void det_remove_anti_fast (MTYPE *storage, double *ts, double *te, int nt, int remid,
		      Matrix * m, Matrix * temp1, Matrix *temp2, Par par);


void det_modify_fast (MTYPE *storage, double *ts, double *te, int nt, int modid,
		      double news, double newe, Matrix * m, Matrix * temp1, Matrix *temp2, Par par);

*/
void propose_time (SEEDTYPE *seed,double *result, double left, double right, double beta);

double myrand(SEEDTYPE *myseed);

double overlap(double as,double ae,double bs,double be);
double overlap_total ( double *ts, double *te, int nt, double start,
		       double end, double beta);
double length_total ( double *ts, double *te, int nt, double beta);

void print_list(int n, int * o_list, double *t_list);

void verify_matrix( int i_o, Mempool * mem, Par par);
