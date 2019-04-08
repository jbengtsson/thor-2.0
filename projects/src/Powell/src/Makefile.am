bin_PROGRAMS  = main main2

main_SOURCES  = main.f newuoa.f
main2_SOURCES = main2.cc newuoa.f

main2_LDADD = -L/usr/lib/x86_64-linux-gnu \
	      -lgfortran

FFLAGS = -g -O2 -Wall -fbounds-check -fdefault-integer-8
