#SUBDIRS      = Powell

bin_PROGRAMS = thor H_min max_iv_match max_iv_uc max_iv_dnu \
	       vlad vlad_2 vlad_3 vlad_4 vlad_5 sls_2 diamond_ii delta \
               delta_lo delta_lo2 quad2d quad3d tune_conf tune_conf_2 \
	       tune_conf_3 tbs get_dnu mpole_scheme poincare_map \
               chk_drv_terms sxt_resp_mat max_4u
#               sext_resp_mat

thor_SOURCES          = thor.cc
H_min_SOURCES         = H_min.cc
max_iv_match_SOURCES  = max_iv_match.cc
max_iv_uc_SOURCES     = max_iv_uc.cc
max_iv_dnu_SOURCES    = max_iv_dnu.cc
vlad_SOURCES          = vlad.cc
vlad_2_SOURCES        = vlad_2.cc
vlad_3_SOURCES        = vlad_3.cc
vlad_4_SOURCES        = vlad_4.cc
vlad_5_SOURCES        = vlad_5.cc
sls_2_SOURCES         = sls_2.cc
diamond_ii_SOURCES    = diamond_ii.cc
delta_SOURCES         = delta.cc
delta_lo_SOURCES      = delta_lo.cc
delta_lo2_SOURCES     = delta_lo2.cc
quad2d_SOURCES        = quad2d.cc
quad3d_SOURCES        = quad3d.cc
tune_conf_SOURCES     = tune_conf.cc
tune_conf_2_SOURCES   = tune_conf_2.cc
tune_conf_3_SOURCES   = tune_conf_3.cc
tbs_SOURCES           = tbs.cc
get_dnu_SOURCES       = get_dnu.cc
mpole_scheme_SOURCES  = mpole_scheme.cc param_type.cc
poincare_map_SOURCES  = poincare_map.cc
chk_drv_terms_SOURCES = chk_drv_terms.cc
sxt_resp_mat_SOURCES  = sxt_resp_mat.cc
max_4u_SOURCES        = max_4u.cc param_type.cc
# sext_resp_mat_SOURCES = sext_resp_mat.c

GSL_INCLUDE = /usr/local/Cellar/gsl/2.8
# GSL_INCLUDE = /net/home/lli/OPAL/OPAL-2.0.1

LIBS	= -L$(THOR_LIB)/thor/lib -lthor \
	  -L$(NUM_REC)/lib -lnum_rec \
	  -lstdc++ \
          -lm \
	  -L$(GSL_INCLUDE)/lib -lgsl \
	  -L/usr/local/opt/gfortran/lib/gcc/14 \
          -lgfortran

#INCLUDES = -I$(THOR_LIB)/thor/inc -I$(NUM_REC)/inc -I$(HOME)/R-2.8.1/include \
#	   -I/usr/local/include -I/usr/include/openmpi

AM_CPPFLAGS = -I$(THOR_LIB)/thor/inc -I$(NUM_REC)/inc  \
              -I$(HOME)/R-2.8.1/include \
	      -I/usr/local/include -I/usr/include/openmpi \
              -I$(GSL_INCLUDE)/include

AM_FFLAGS   = -g -O0 -Wall -fbounds-check

AM_CXXFLAGS = -g -O2 -Wall -fPIC -std=c++11

UNLIMIT_RESOURCES = datasize memorysize stacksize
