lib_LIBRARIES = libthor.a

libthor_a_SOURCES = thor_lib.cc \
                    ../../TPSA/TPSALib.f ../../TPSA/LieLib.f

#INCLUDES = -I../inc -I$(NUM_REC)/inc
GSL_INCLUDE = -I/usr/local/Cellar/gsl/2.7.1/include
# GSL_INCLUDE = -I/net/home/lli/OPAL/OPAL-2.0.1/include
AM_CPPFLAGS = $(GSL_INCLUDE) \
	      -I../inc -I$(NUM_REC)/inc
#LIBS     = -L/home/bengtsson/gcc-4.2/lib

FFLAGS      = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=medium \
	      -w -fallow-argument-mismatch
#AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
