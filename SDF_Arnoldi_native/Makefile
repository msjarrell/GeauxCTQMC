#MKLROOT = /opt/intel/mkl
MKLROOT=/usr/local/compilers/Intel/cluster_studio_xe_2013.1.046/composer_xe_2013_sp1.2.144/mkl

OPTFLAGS = -O3
xOPTFLAGS = 

TIMING ?=1
FOR_MIC ?= 0
USE_GCC ?= 0

ifeq ($(USE_GCC),1)
CXX = g++
else
CXX = icpc
endif

SRC_FILES = main.cpp mat_exp_krylov.cpp Arnoldi.cpp common.cpp \
	SDF_SPMV_dc_func_n.cpp SDF_SPMV_dc_func_p.cpp

ISRC_FILES = initialize_I.cpp SDF_SPMV_init_I.cpp SDF_SPMV_init_txt_I.cpp \
	CSR2DIA_I.cpp DIA2SDF_I.cpp SDF_SPMV_I.cpp

LIB_SRC = mat_exp_krylov.cpp Arnoldi.cpp common.cpp \
	SDF_SPMV_dc_func_n.cpp SDF_SPMV_dc_func_p.cpp

ILIB_SRC = SDF_SPMV_init_I.cpp SDF_SPMV_init_txt_I.cpp \
	CSR2DIA_I.cpp DIA2SDF_I.cpp SDF_SPMV_I.cpp

ALL_SRC_FILES = $(SRC_FILES) $(ISRC_FILES)

OBJ_FILES = $(SRC_FILES:.cpp=.o)
IOBJ_FILES = $(ISRC_FILES:.cpp=.o)
# $(INCLUDE)

ifeq ($(USE_GCC),1)
CXXFLAGS = $(OPTFLAGS) -g  -Wall -Wno-parentheses
else
CXXFLAGS = $(OPTFLAGS) -g  -DHAVE_MKL
endif

ifeq ($(TIMING),1)
CXXFLAGS += -DTIMING
endif

ifeq ($(USE_GCC),1)
ASM_GEN_FLAGS = --verbose-asm
else
ASM_GEN_FLAGS = -fverbose-asm
endif

ifeq ($(USE_GCC),1)
LIBS += -lgfortran -llapack -lblas -lrt
else
ifeq ($(FOR_MIC),1)
LIBS += -lrt -L/opt/intel/lib/mic -I$(MKLROOT)/include \
              -Wl,--start-group $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
              $(MKLROOT)/lib/mic/libmkl_sequential.a \
              $(MKLROOT)/lib/mic/libmkl_core.a -Wl,--end-group
else
LIBS += -lrt -lmkl_intel_lp64 -lmkl_core -lmkl_sequential 
endif
endif

ICXXFLAGS := $(CXXFLAGS)

ifeq ($(FOR_MIC),1)
CXXFLAGS += -mmic
endif

ILINKFLAGS = $(ICXXFLAGS)

LINKFLAGS = $(CXXFLAGS)

DEFAULT: sdfexpn_I.so main sdfexpn.so


# Include the dependencies that were created by %.d rule.
#
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_SRC_FILES:.cpp=.d)
endif
#

# Prepare file holding dependencies, to be included in this file.
#
%.d: %.cpp Makefile
	@set -e; rm -f $@; \
	g++ -DMAKE_DEPEND -M $(ICXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%_I.o: %_I.cpp Makefile
	$(CXX) $(ICXXFLAGS) -S $< $(ASM_GEN_FLAGS)
	$(CXX) $(ICXXFLAGS) -c $<

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -S $< $(ASM_GEN_FLAGS)
	$(CXX) $(CXXFLAGS) -c $<


initialize: $(IOBJ_FILES)
	$(CXX) $(ILINKFLAGS) -o $@ $(IOBJ_FILES)

SDF_SPMV_dc_func_n.cpp: initialize
	./initialize

SDF_SPMV_dc_func_p.cpp: initialize
	./initialize

main: $(OBJ_FILES)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJ_FILES) $(LIBS)

sdfexpn.so: $(LIB_SRC) Makefile
	$(CXX) $(LINKFLAGS) -shared -fPIC -o $@ $(LIB_SRC) $(LIBS)

sdfexpn_I.so: $(ILIB_SRC) Makefile
	$(CXX) $(ILINKFLAGS) -shared -fPIC -o $@ $(ILIB_SRC)

clean:
	rm -f $(IOBJ_FILES) $(OBJ_FILES) *.o *.s *.d \
	main initialize SDF_SPMV_dc_func_p.cpp SDF_SPMV_dc_func_n.cpp vecd.txt sdfexpn.so sdfexpn_I.so
