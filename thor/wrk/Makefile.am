bin_PROGRAMS =  vlad

vlad_SOURCES = vlad.cc

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
