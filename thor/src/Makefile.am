lib_LIBRARIES = libthor.a

libthor_a_SOURCES = thor_lib.cc \
                    ../../TPSA/TPSALib.f ../../TPSA/LieLib.f

#INCLUDES = -I../inc -I$(NUM_REC)/inc
AM_CPPFLAGS = -I../inc -I$(NUM_REC)/inc
#LIBS     = -L/home/bengtsson/gcc-4.2/lib

#FFLAGS   = -g -O2 -Wall -fbounds-check -i4 -mcmodel=medium
#FFLAGS   = -g -O2 -Wall -fbounds-check -mcmodel=large -fstack-arrays
#FFLAGS   = -g -O2 -Wall -fbounds-check -mcmodel=large -frecursive
FFLAGS   = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large
#AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend
AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC
