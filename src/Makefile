.SUFFIXES: .c .h

#CC = /usr/local/compilers/Intel/parallel_studio_xe_2015/composer_xe_2015.0.090/bin/intel64/icc 
#CC = /usr/local/compilers/Intel/parallel_studio_xe_2015/composer_xe_2015.0.090/bin/intel64/icpc
#CC = /usr/local/compilers/Intel/cluster_studio_xe_2013.1.046/composer_xe_2013_sp1.2.144/bin/intel64/icpc
CC= ~/intel/composer_xe_2015.3.187/bin/intel64/icpc

MPICC= /usr/local/compilers/Intel/cluster_studio_xe_2013.1.046/impi/4.1.3.048/intel64/bin/mpiicc
#INC = -I/usr/local/compilers/Intel/composer_xe_2013.5.192/mkl/include/
INC= -I/usr/local/compilers/Intel/parallel_studio_xe_2015/composer_xe_2015.0.090/mkl/include/
#LIB_MIC = -L/usr/local/compilers/Intel/composer_xe_2013.5.192/mkl/lib/mic/ #-L/usr/local/packages/lapack/3.4.0/gcc-4.4.6/lib
LIB_MIC= -L/usr/local/compilers/Intel/parallel_studio_xe_2015/composer_xe_2015.0.090/mkl/lib/mic/
#LIB_CPU = -L/usr/local/compilers/Intel/composer_xe_2013.5.192/mkl/lib/intel64/ #-L/usr/local/packages/lapack/3.4.0/gcc-4.4.6/lib
#LIB_CPU= -L/usr/local/compilers/Intel/parallel_studio_xe_2015/composer_xe_2015.0.090/mkl/lib/intel64/
LIB_CPU= -L~/intel/composerxe/mkl/lib/intel64/

#CC= /opt/intel/composer_xe_2013_sp1.1.106/bin/intel64/icc
#INC = -I/opt/intel/composer_xe_2013_sp1.1.106/mkl/include/
#LIB_MIC = -L/opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/mic/ #-L/usr/local/packages/lapack/3.4.0/gcc-4.4.6/lib
#LIB_CPU = -L/opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64/ #-L/usr/local/packages/lapack/3.4.0/gcc-4.4.6/lib

#LIB = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm  #-llapack -lgfortran -lblas
LIB = -mkl=parallel -lpthread -lm  #-llapack -lgfortran -lblas


OPTFLAG2 = -O2 #-vec-report3
OPT_FLAG = -O2

CPU_FLAGS = $(INC) $(LIB) $(LIB_CPU) $(ARCHFLAG) 
MIC_FLAGS = $(INC) $(LIB) $(LIB_MIC) $(ARCHFLAG)

DEBUG_FLAG = #-g #-DDEBUG_SHIFT
VERBOSE= #-DDEBUG_VERBOSE  -DDEBUG_VERBOSE_SHIFT 

LDFLAGS =  -L./SDF_exp_API_native/ ./SDF_exp_API_native/sdfexpn.so

src_det =  main_matrix.c init.c calc_green_func.c
exec_det= calc_det

#mpi_src =  accept.c  proposal.c weight.c calc_green_func.c main.c utility.c measurement.c mpi_wrapper.c det.c
#mpi_src =  accept.c  proposal.c weight_dhm.c calc_green_func.c main.c utility.c measurement.c mpi_wrapper.c det.c
mpi_src =  accept.c  proposal.c weight_aim.c calc_green_func.c main.c utility.c measurement.c mpi_wrapper.c det.c
#st_src =  accept.c  proposal.c weight.c calc_green_func.c main.c utility.c measurement.c single_wrapper.c det.c 
#st_src =  accept.c  proposal.c weight_dhm.c calc_green_func.c main.c utility.c measurement.c single_wrapper.c det.c 
st_matrix_src =  accept.c  proposal.c weight_aim.c calc_green_func.c main.c utility.c measurement.c single_wrapper.c det.c 
st_src =  accept.c  proposal.c weight.c calc_green_func.c main.c utility.c measurement.c single_wrapper.c det.c 
exec = ctqmc_cpu
exec_mic = ctqmc_mic
exec_g = ctqmc_g
exec_g_matrix = ctqmc_g_matrix
exec_g_mic = ctqmc_g_mic
exec_g_sdf = ctqmc_g_sdf

test_src = hamiltonian_test.c
exec_test = hamiltonian_test
exec_dmft = ctqmc_dmft

dhm_src =  accept.c  proposal.c weight_dhm.c calc_green_func.c main.c utility.c measurement.c single_wrapper.c det.c 

exe_dhm = ctqmc_dhm

#SDF_SRC_FILES = main.cpp mat_exp_krylov.cpp Arnoldi.cpp common.cpp \
	SDF_SPMV_dc_func_n.cpp SDF_SPMV_dc_func_p.cpp
sdf_f=./SDF_Arnoldi_native
#SDF_O = $(sdf_f)/mat_exp_krylov.o $(sdf_f)/Arnoldi.o $(sdf_f)/common.o \
	$(sdf_f)/SDF_SPMV_dc_func_n.o $(sdf_f)/SDF_SPMV_dc_func_p.o


default:  $(exe_dmft)
#$(exec_g) $(exec_g_matrix)
#$(exec) $(exec_mic) SDF

SDF: 
	make -C SDF_exp_API_native

$(exe_dhm): $(dhm_src)
	$(CC) $(OPT_FLAG)$(DEBUG_FLAG) $(VERBOSE) -DDMFT_LOOP -DUSE_MATRIX -DHAVE_MKL  -fopenmp  $(dhm_src)  $(CPU_FLAGS) -o $(exe_dhm)
#	$(CC) $(OPT_FLAG)$(DEBUG_FLAG) $(VERBOSE) -DUSE_MATRIX -DHAVE_MKL  -fopenmp  $(dhm_src)  $(CPU_FLAGS) -o $(exe_dhm)



$(exec_g): $(st_src)
#	$(CC) -DSDF -DHAVE_MKL -g  -fopenmp  $(st_src)  $(CPU_FLAGS) -o $(exec_g) $(LDFLAGS) -lsvml -limf -lirng -lintlc -I $(sdf_f)
	$(CC) $(DEBUG_FLAG) $(VERBOSE) -DMEASURE_LEG -DDMFT_LOOP -DHAVE_MKL -g  -fopenmp  $(st_src)  $(CPU_FLAGS) -o $(exec_g)
# $(LDFLAGS) -lsvml -limf -lirng -lintlc -I $(sdf_f)
#	$(CC) -DSDF -DHAVE_MKL -g  -fopenmp  $(st_src)  $(CPU_FLAGS) -o $(exec_g) $(SDF_O) $(LDFLAGS) -lsvml -limf -lirng -lintlc -I ./SDF_Arnoldi_native
#	$(CC) -DHAVE_MKL -g  -fopenmp  $(st_src)  $(CPU_FLAGS) -o $(exec_g) $(SDF_O) $(LDFLAGS) -lsvml -limf -lirng -lintlc -I ./SDF_Arnoldi_native

$(exec_g_matrix): $(st_src_matrix)
#	$(CC) -DSDF -DHAVE_MKL -g  -fopenmp  $(st_src)  $(CPU_FLAGS) -o $(exec_g) $(LDFLAGS) -lsvml -limf -lirng -lintlc -I $(sdf_f)
#	$(CC)  $(DEBUG_FLAG) $(VERBOSE) -DDMFT_LOOP -DUSE_MATRIX -DHAVE_MKL -g  -fopenmp  $(st_matrix_src)  $(CPU_FLAGS) -o $(exec_g_matrix) 
	$(CC)  $(DEBUG_FLAG) $(VERBOSE) -DMEASURE_GIO -DMEASURE_LEG -DUSE_MATRIX -DHAVE_MKL -g  -fopenmp  $(st_matrix_src)  $(CPU_FLAGS) -o $(exec_g_matrix) 
#	$(CC)  $(DEBUG_FLAG) $(VERBOSE) -DUSE_MATRIX -DHAVE_MKL -g  -fopenmp  $(st_matrix_src)  $(CPU_FLAGS) -o $(exec_g_matrix) 
#$(LDFLAGS) -lsvml -limf -lirng -lintlc -I $(sdf_f)



#$(exec_g_sdf): $(st_src_matrix)
#	$(CC)  $(DEBUG_FLAG) $(VERBOSE) -DSDF -DMEASURE_GIO -DMEASURE_LEG -DUSE_MATRIX -DHAVE_MKL -g  -fopenmp  $(st_src_matrix)  $(CPU_FLAGS) -o $(exec_g_matrix) $(LDFLAGS) -lsvml -limf -lirng -lintlc -I $(sdf_f)



$(exec_dmft): $(st_src_matrix)
	$(CC)  $(DEBUG_FLAG) $(VERBOSE) -DMEASURE_LEG -DDMFT_LOOP -DUSE_MATRIX -DHAVE_MKL -g  -fopenmp  $(st_matrix_src)  $(CPU_FLAGS) -o $(exec_dmft) 




$(exec_g_sdf): $(st_src_matrix) SDF
	$(CC)  $(DEBUG_FLAG) $(VERBOSE) -DSDF -DUSE_MATRIX -DHAVE_MKL -g  -fopenmp  $(st_matrix_src)  $(CPU_FLAGS) -o $(exec_g_sdf) $(SDF_O) $(LDFLAGS) -lsvml -limf -lirng -lintlc -I ./SDF_Arnoldi_native


$(exec_test): $(test_src)
	$(CC) -g  -fopenmp $(test_src)  $(CPU_FLAGS) -o $(exec_test)


$(exec_g_mic): $(st_src)
#	$(CC) -g -mmic -fopenmp $(st_src)  $(MIC_FLAGS) -o $(exec_g_mic)
	$(CC) -DSDF -DHAVE_MKL -g -mmic  -fopenmp  $(st_src)  $(MIC_FLAGS) -o $(exec_g_mic) $(SDF_O) $(LDFLAGS) -lsvml -limf -lirng -lintlc -I ./SDF_Arnoldi_native -lrt

$(exec_mic):
	$(MPICC) -fopenmp -mmic $(mpi_src) $(OPTFLAG2) $(MIC_FLAGS) -o $(exec_mic)

$(exec): 
	$(MPICC) -fopenmp $(mpi_src) $(OPTFLAG2) $(CPU_FLAGS) -o $(exec)

$(exec_det):
	$(CC) $(src_det) $(OPTFLAG2) $(CFLAGS) -o $(exec_det)


clean:
	rm -rf *.o $(exec_g)			
	make -C ./SDF_Arnoldi_native/ clean
	make -C ./SDF_exp_API_native/ clean
