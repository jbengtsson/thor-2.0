bin_PROGRAMS      = thor h_rnd optics tst_field get_SRM H_min H_zero nl_opt \
                    fel tst_cg nl_cg nl_pwl nl_ns ctrl_a3
#                    nl_opt_mpi

thor_SOURCES      = thor.cc
h_rnd_SOURCES     = h_rnd.cc
optics_SOURCES    = optics.cc
tst_field_SOURCES = tst_field.cc
get_SRM_SOURCES   = get_SRM.cc
H_min_SOURCES     = H_min.cc
H_zero_SOURCES    = H_zero.cc
nl_opt            = nl_opt.cc
fel               = fel.cc
tst_cg            = tst_cg.cc
nl_cg             = nl_cg.cc
nl_pwl            = nl_pwl.cc
nl_sd             = nl_sd.cc
ctrl_a3_SOURCES   = ctrl_a3.cc
nl_opt_mpi        = nl_opt_mpi.cc

#nsga2_SOURCES     = nsga2.cc

#fit_cell_SOURCES  = fit_cell.cc

LIBS	= -L$(THOR_LIB)/thor/lib -lthor \
	  -L$(NUM_REC)/lib -lnum_rec \
	  -L/usr/lib \
	  -lstdc++ -lgfortran \
          -lm \
	  -lgsl -lgslcblas
#	  -lmpich \
#	  -lmpl
#          -lg2c 
#	  -lstdc++ -ltermcap -ldl \
#	  -L$(HOME)/R-2.8.1/src/main -L$(HOME)/R-2.8.1/lib
#          -lR -lRblas -lRlapack
#	  -L$(HOME)/readline-5.0
#	  -lreadline -lhistory

INCLUDES = -I$(THOR_LIB)/thor/inc -I$(NUM_REC)/inc -I$(HOME)/R-2.8.1/include \
	   -I/usr/local/include -I/usr/include/openmpi

CXXFLAGS = -g -O2 -Wall -Wno-non-template-friend
#CXXFLAGS = -g -Wl,-stack_size -Wl,0x40000000 -O2 -Wall -Wno-non-template-friend
UNLIMIT_RESOURCES = datasize memorysize stacksize
