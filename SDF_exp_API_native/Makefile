ICCROOT = /opt/intel/composer_xe_2013_sp1.1.106/compiler
#MKLROOT = /opt/intel/mkl
MKLROOT=/usr/local/compilers/Intel/cluster_studio_xe_2013.1.046/composer_xe_2013_sp1.2.144/mkl

SDFEXP_DIR = ../SDF_Arnoldi_native
SDFEXP_INC = $(SDFEXP_DIR)
SDFEXP_LIB = $(SDFEXP_DIR)

CXXFLAGS = -Wall -I $(SDFEXP_INC)
ILDFLAGS =  -L$(SDFEXP_LIB) $(SDFEXP_LIB)/sdfexpn_I.so

LDFLAGS =  -L./ ./sdfexpn.so

LINKFLAGS = $(CXXFLAGS) -O3 -g  -DHAVE_MKL #-mmic

LIB_SRC = SDF_SPMV_dc_func_n.cpp SDF_SPMV_dc_func_p.cpp \
	$(SDFEXP_LIB)/mat_exp_krylov.cpp $(SDFEXP_LIB)/Arnoldi.cpp \
	$(SDFEXP_LIB)/common.cpp 

#LIBS += -lrt -L/opt/intel/lib/mic -I$(MKLROOT)/include \
	-Wl,--start-group $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
        $(MKLROOT)/lib/mic/libmkl_sequential.a \
        $(MKLROOT)/lib/mic/libmkl_core.a -Wl,--end-group 

LIBS += -lrt -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

CXX = /usr/local/compilers/Intel/cluster_studio_xe_2013.1.046/composer_xe_2013_sp1.2.144/bin/intel64/icpc -DHAVE_MKL

DEFAULT: main

$(SDFEXP_INC)/sdfexpn_I.so:
	$(MAKE) -C $(SDFEXP_DIR)

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -c $<

SDF_SPMV_dc_func_n.cpp: main_I
	./main_I

SDF_SPMV_dc_func_p.cpp: main_I
	./main_I

main_I: main_I.o $(SDFEXP_INC)/sdfexpn_I.so
	$(CXX) $< -o $@ $(CXXFLAGS) $(ILDFLAGS)

sdfexpn.so: $(LIB_SRC) Makefile
	$(CXX) $(LINKFLAGS) -shared -fPIC -o $@ $(LIB_SRC) $(LIBS)

main: $(LIB_SRC) Makefile sdfexpn.so
	$(CXX) $(LINKFLAGS) -o $@ main.cpp $(LDFLAGS) -lsvml -limf -lirng -lintlc

clean:
	rm -f *.o *.s *.d *.so SDF_SPMV_dc_func_n.cpp SDF_SPMV_dc_func_p.cpp main_I main vecd.txt 