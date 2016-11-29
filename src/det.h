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


void det_shift_fast (MTYPE *storage, double *ts, double *te, int nt, int modid,
		      double news, double newe, Matrix * m, Matrix * temp1, Matrix *temp2, Par par);


void construct_matrix(double *ts, double *te, int nt, Matrix *temp, Par par);
