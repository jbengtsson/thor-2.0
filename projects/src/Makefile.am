bin_PROGRAMS = H_min max_iv_match max_iv_uc max_iv_dnu \
	       vlad vlad_2 vlad_3 vlad_4 vlad_5 sls_2 diamond_ii delta \
               delta_lo delta_lo2 quad2d quad3d tune_conf tune_conf_2 \
	       tune_conf_3 tbs

H_min_SOURCES        = H_min.cc
max_iv_match_SOURCES = max_iv_match.cc
max_iv_uc_SOURCES    = max_iv_uc.cc
max_iv_dnu_SOURCES   = max_iv_dnu.cc
vlad_SOURCES         = vlad.cc
vlad_2_SOURCES       = vlad_2.cc
vlad_3_SOURCES       = vlad_3.cc
vlad_4_SOURCES       = vlad_4.cc
vlad_5_SOURCES       = vlad_5.cc
sls_2_SOURCES        = sls_2.cc
diamond_ii_SOURCES   = diamond_ii.cc
delta_SOURCES        = delta.cc
delta_lo_SOURCES     = delta_lo.cc
delta_lo2_SOURCES    = delta_lo2.cc
quad2d_SOURCES       = quad2d.cc
quad3d_SOURCES       = quad3d.cc
tune_conf_SOURCES    = tune_conf.cc
tune_conf_2_SOURCES  = tune_conf_2.cc
tune_conf_3_SOURCES  = tune_conf_3.cc
tbs_SOURCES          = tbs.cc

LIBS	= -L$(THOR_LIB)/thor/lib -lthor \
	  -L$(NUM_REC)/lib -lnum_rec \
	  -L/usr/lib/x86_64-linux-gnu \
	  -lstdc++ -lgfortran \
          -lm
#	  -lgsl -lgslcblas
#	  -lmpich \
#	  -lmpl
#          -lg2c 
#	  -lstdc++ -ltermcap -ldl \
#	  -L$(HOME)/R-2.8.1/src/main -L$(HOME)/R-2.8.1/lib
#          -lR -lRblas -lRlapack
#	  -L$(HOME)/readline-5.0
#	  -lreadline -lhistory

#INCLUDES = -I$(THOR_LIB)/thor/inc -I$(NUM_REC)/inc -I$(HOME)/R-2.8.1/include \
#	   -I/usr/local/include -I/usr/include/openmpi
AM_CPPFLAGS = -I$(THOR_LIB)/thor/inc -I$(NUM_REC)/inc  \
              -I$(HOME)/R-2.8.1/include \
	      -I/usr/local/include -I/usr/include/openmpi

#AM_CXXFLAGS = -g -O2 -Wall -Wno-non-template-friend -fPIC -Wl,--no-relax
AM_CXXFLAGS = -g -O2 -Wall -Wno-non-template-friend -fPIC
#CXXFLAGS = -g -Wl,-stack_size -Wl,0x40000000 -O2 -Wall -Wno-non-template-friend
UNLIMIT_RESOURCES = datasize memorysize stacksize
