#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = pddrive Bogus


Bogus_SOURCES = $(srcdir)/Bogus

Bogus$(EXEEXT): $(srcdir)/Bogus
	cp  $(srcdir)/Bogus Bogus$(EXEEXT) 
	cp  $(srcdir)/g20.rua.in g20.rua 

pddrive_CFLAGS = -I$(top_srcdir)/SRC

pddrive_SOURCES = $(srcdir)/pddrive.c $(srcdir)/dcreate_matrix.c

pddrive_DEPENDENCIES = \
	$(top_builddir)/SRC/libsuperludist.a 

pddrive_LDADD = $(top_builddir)/SRC/libsuperludist.a \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS) 

