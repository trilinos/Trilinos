#@HEADER
#@HEADER

EXEEXT = .exe
# Nothing should be done for an mpi build
if HAVE_MPI
else

noinst_PROGRAMS = epetraTestAll

epetraTestAll_SOURCES = $(srcdir)/epetraTestAll

epetraTestAll$(EXEEXT): $(srcdir)/epetraTestAll
	cp $(srcdir)/epetraTestAll epetraTestAll$(EXEEXT) 

endif
