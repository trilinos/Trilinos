
srcdir = .
top_srcdir = ..
top_srcdir = /home/dpdiach/Trilinos/packages/nox

INCLUDE = -I$(top_srcdir)/src -I$(top_srcdir)/src-lapack -I/. \
	  -I../../../test/utils

LIBS	= $(top_srcdir)/src-lapack/libnoxlapack.a \
	  $(top_srcdir)/src/libnox.a -llapack -lblas \
	  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2 \
	  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/../../.. \
	  -lfrtbegin -lg2c -lm -lgcc_s

FLAGS	= -O2 -DHAVE_CONFIG_H -DWITH_PRERELEASE  $(INCLUDE)

OBJECT	= Rosenbrock.o Utils.o

CC	= g++

rosenbrock:	$(OBJECT)
		$(CC) -o rosenbrock $(OBJECT) $(LIBS)

.C.o:				
		$(CC) $(FLAGS) $(INCLUDES) -c $*.C

