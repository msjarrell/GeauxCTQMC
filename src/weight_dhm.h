void weight_add(Mempool *mempool, int up, int newid, double news, double newe, Par par);
void weight_remove(Mempool *mempool, int up, int remid, Par par);
void weight_shift(Mempool *mempool, int up, int modid, double news, double newe, Par par);
void weight_add_anti(Mempool *mempool, int up, int newid, double news, double newe, Par par);
void weight_remove_anti(Mempool *mempool, int up, int remid, Par par);



void exponential_matrix(int m, double *v, double t,double *expvt);

void generate_shift_list(int n,int *op_list,double *time_list, int *tmp_op_list, double * tmp_time_list, int op, double oldt, double newt);

void generate_add_list(int n,int *op_list,double *time_list, int *tmp_op_list, double * tmp_time_list, int i_o,  double news, double newe);

void generate_rem_list(int n,int *op_list,double *time_list, int *tmp_op_list, double * tmp_time_list, int i_o, double rems,double reme);
  
double weight_loc(int n, int *op_list,double *time_list, Par par);

void mmmul(double *m1, double *m2);//multiply matrices m1 and m2, and store result in m2

void copy_matrix( double *m1, double *m2);

void init_hamiltonian(Par *par);
