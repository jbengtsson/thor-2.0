AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libTPSALib.a libLieLib.a

libTPSALib_a_SOURCES = TPSALib.f

libLieLib_a_SOURCES  = LieLib.f

#FFLAGS = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large 
# Option "large" is not supported for the memory model.
FFLAGS = \
	-g -Og -Wall -fbounds-check -fdefault-integer-8 -mcmodel=medium \
	-w -fallow-argument-mismatch

#AM_FFLAGS = -g -O2 -Wall -fbounds-check -fPIC -mcmodel=medium
#AM_FFLAGS = -g -O2 -Wall -fbounds-check -fPIC -mcmodel=large -fstack-arrays
AM_FFLAGS = -g -Og -fbacktrace -Wall -fbounds-check -fPIC
